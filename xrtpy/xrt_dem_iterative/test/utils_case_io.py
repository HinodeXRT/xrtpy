from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import List, Tuple

import numpy as np
import pandas as pd

from scipy.io import readsav

#from xrtpy.response.tools import generate_temperature_responses
from xrtpy.response.tools.multi_filter_response import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative



def case_dir(case_name: str) -> Path:
    return Path(__file__).parent / "data" / "cases" / case_name


@dataclass(frozen=True)
class MonteCarloInputs:
    df: pd.DataFrame
    filters: List[str]
    mc_intensities: np.ndarray  # shape: (n_runs, n_filters)


def read_mc_intensities_csv(csv_path: str | Path) -> MonteCarloInputs:
    """
    Read Monte Carlo intensities CSV exported from IDL.

    Expected format:
        run,<filter1>,<filter2>,...
        1, ...
        2, ...
        ...

    Returns:
        MonteCarloInputs with:
            df: DataFrame (run + filter columns)
            filters: list of filter column names
            mc_intensities: float ndarray of shape (n_runs, n_filters)
    """
    csv_path = Path(csv_path)
    if not csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {csv_path}")

    df = pd.read_csv(csv_path)

    # basic structure checks
    if df.shape[1] < 2:
        raise ValueError(
            f"CSV must have at least 2 columns (run + filters). Got columns: {list(df.columns)}"
        )

    # Ensure we have a 'run' column (IDL sometimes writes 'run' exactly, but be safe)
    if "run" not in df.columns:
        # try case-insensitive match
        matches = [c for c in df.columns if c.strip().lower() == "run"]
        if not matches:
            raise ValueError(f"CSV missing a 'run' column. Found columns: {list(df.columns)}")
        df = df.rename(columns={matches[0]: "run"})

    # normalize run
    df["run"] = pd.to_numeric(df["run"], errors="raise").astype(int)

    # Sort by run just in case, and reset index
    df = df.sort_values("run").reset_index(drop=True)

    # Require runs to be 1..N (strict on purpose for testing)
    expected = np.arange(1, len(df) + 1)
    if not np.array_equal(df["run"].to_numpy(), expected):
        raise ValueError(
            f"Run column must be consecutive 1..N.\n"
            f"Expected: {expected[:10]}...\n"
            f"Got:      {df['run'].to_numpy()[:10]}..."
        )

    # filter columns
    filters = [c for c in df.columns if c != "run"]
    if len(filters) == 0:
        raise ValueError("No filter columns found (only 'run' present).")

    # Convert all filter columns to float
    for c in filters:
        df[c] = pd.to_numeric(df[c], errors="raise").astype(float)

    mc_intensities = df[filters].to_numpy(dtype=float) #(n_runs, n_filters)

    # No NaNs - but shoudn't have any
    if not np.isfinite(mc_intensities).all():
        bad = np.where(~np.isfinite(mc_intensities))
        raise ValueError(f"Found non-finite intensity values at indices: {bad}")

    return MonteCarloInputs(df=df, filters=filters, mc_intensities=mc_intensities)



@dataclass(frozen=True)
class DemBatchResult:
    filters: list[str]
    logT: np.ndarray# (nT,)
    dem_runs: np.ndarray# (n_runs, nT)
    modeled_runs: np.ndarray # (n_runs, n_filters)
    chisq_runs: np.ndarray# (n_runs,)


def run_dem_for_mc_csv(
    csv_path: str | Path,
    observation_date: str,
    *,
    intensity_errors: np.ndarray | None = None,
    minimum_bound_temperature: float = 5.5,
    maximum_bound_temperature: float = 8.0,
    logarithmic_temperature_step_size: float = 0.1,
) -> DemBatchResult:
    """
    Run XRTpy DEM for each row in an IDL-exported MC intensity CSV.

    Notes:
    - We run monte_carlo_runs=0 because the CSV already provides the perturbed intensities.
    - This makes the Python/IDL comparison deterministic.
    """
    mc = read_mc_intensities_csv(csv_path)
    filters = mc.filters
    mc_intensities = mc.mc_intensities  # (n_runs, n_filters)

    # Generate temperature responses once
    responses = generate_temperature_responses(filters, observation_date)

    # Create a solver "template" once (we'll reuse and only replace intensities)
    solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=mc_intensities[0],
        temperature_responses=responses,
        intensity_errors=intensity_errors,  # optional; can be None - Make sure current xrtpy-dem has it set- NOTEFORJOY
        monte_carlo_runs=0,
        minimum_bound_temperature=minimum_bound_temperature,
        maximum_bound_temperature=maximum_bound_temperature,
        logarithmic_temperature_step_size=logarithmic_temperature_step_size,
    )

    # Build grid/response matrix once
    solver.create_logT_grid()
    solver._interpolate_responses_to_grid()

    n_runs, n_filters = mc_intensities.shape
    nT = len(solver.logT)

    dem_runs = np.zeros((n_runs, nT), dtype=float)
    modeled_runs = np.zeros((n_runs, n_filters), dtype=float)
    chisq_runs = np.zeros((n_runs,), dtype=float)

    # Run each case deterministically
    for i in range(n_runs):
        dem_i, modeled_i, chi2_i, _ = solver._solve_single_dem(
            observed_intensities_vals=mc_intensities[i]
        )
        dem_runs[i] = dem_i
        modeled_runs[i] = modeled_i
        chisq_runs[i] = chi2_i

    return DemBatchResult(
        filters=filters,
        logT=np.array(solver.logT, dtype=float),
        dem_runs=dem_runs,
        modeled_runs=modeled_runs,
        chisq_runs=chisq_runs,
    )


@dataclass(frozen=True)
class IDLDemResult:
    logT: np.ndarray# (nT,)
    dem_runs: np.ndarray# (n_runs, nT)   (runs exclude base or include base depending on file)
    dem_base: np.ndarray # (nT,)
    n_runs: int


def _to_1d(x) -> np.ndarray:
    x = np.array(x)
    return x.ravel()


def _ensure_runs_by_T(dem: np.ndarray, logT: np.ndarray) -> np.ndarray:
    """
    Return DEM shaped (n_runs, nT).

    Handles common IDL orientations:
        - (nT, nRuns+1)
        - (nRuns+1, nT)
        - sometimes (nT, nRuns) etc.

    We infer nT from logT length.
    """
    dem = np.array(dem)
    nT = len(logT)

    if dem.ndim != 2:
        raise ValueError(f"Expected 2D DEM array from IDL, got shape {dem.shape}")

    a, b = dem.shape

    # Case 1: columns are runs, rows are temperature
    if a == nT:
        # dem is (nT, nRunsX) -> transpose to (nRunsX, nT)
        return dem.T

    # Case 2: rows are runs, columns are temperature
    if b == nT:
        # dem already (nRunsX, nT)
        return dem

    raise ValueError(
        f"Cannot infer DEM orientation: dem shape={dem.shape}, logT length={nT}"
    )


def load_idl_dem_sav(sav_path: str | Path) -> IDLDemResult:
    """
    Load IDL .sav file produced by your IDL DEM script and return:
        - logT (nT,)
        - dem_base (nT,)
        - dem_runs (n_runs, nT)  (MC only, base removed)
    """
    
    data = readsav(str(sav_path), python_dict=True)

    # Be flexible about key names
    # Your newer IDL saver used logT_out/dem_out; older might be logt/dem_mc/dem
    logT_key_candidates = ["logT_out", "logt_out", "logt", "logT", "logT_mc", "logT_idl"]
    dem_key_candidates = ["dem_out", "dem_mc", "dem", "dem0"]

    logT = None
    for k in logT_key_candidates:
        if k in data:
            logT = _to_1d(data[k])
            break
    if logT is None:
        raise KeyError(f"Could not find logT in .sav. Keys: {sorted(data.keys())}")

    dem = None
    for k in dem_key_candidates:
        if k in data:
            dem = np.array(data[k])
            break
    if dem is None:
        raise KeyError(f"Could not find DEM in .sav. Keys: {sorted(data.keys())}")

    dem_runs_all = _ensure_runs_by_T(dem, logT)  # (nRunsX, nT)

    # Convention: in IDL output, column/run 0 is usually the base DEM
    dem_base = dem_runs_all[0].copy()
    dem_mc = dem_runs_all[1:].copy()

    return IDLDemResult(
        logT=logT,
        dem_base=dem_base,
        dem_runs=dem_mc,
        n_runs=dem_mc.shape[0],
    )

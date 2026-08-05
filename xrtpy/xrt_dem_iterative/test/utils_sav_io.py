"""
utils_sav_io.py
===============
Shared utilities for loading IDL DEM .sav files and parsing their filenames.

Filename convention expected:
    xrt_IDL_dem_<DATE>_<Filter1><Intensity1>_<Filter2><Intensity2>_..._.sav

Example:
    xrt_IDL_dem_20071213T0401_Bemed603.875886_Bethin150.921435_Alpoly2412.34_.sav

Compact filter name  XRTpy filter name mapping is handled automatically.
Will be updated once I start testing with MC.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.io import readsav

# ---------------------------------------------------------------------------
# Filter name mapping
# ---------------------------------------------------------------------------
_FILTER_MAP: dict[str, str] = {
    "AlpolyAlmesh": "Al-poly/Al-mesh",
    "AlpolyTipoly": "Al-poly/Ti-poly",
    "AlpolyAlthick": "Al-poly/Al-thick",
    "AlpolyBethick": "Al-poly/Be-thick",
    "CpolyTipoly": "C-poly/Ti-poly",
    "Bemed": "Be-med",
    "Bethin": "Be-thin",
    "Bethick": "Be-thick",
    "Alpoly": "Al-poly",
    "Tipoly": "Ti-poly",
    "Almesh": "Al-mesh",
    "Almed": "Al-med",
    "Althick": "Al-thick",
    "Cpoly": "C-poly",
}
_FILTER_KEYS_SORTED = sorted(_FILTER_MAP, key=len, reverse=True)


@dataclass(frozen=True)
class SavCase:
    """
    All information needed to run one IDL vs XRTpy comparison.

    Attributes
    ----------
    sav_path : Path
        Path to the IDL .sav file.
    observation_date : str
        ISO-8601 observation date, e.g. ``"2008-01-04T11:04:26"``.
    filters : list[str]
        XRTpy filter names, in the same order as intensities.
    intensities : list[float]
        Observed intensities in DN/s/pix, matching filters.
    label : str
        Short human-readable label derived from the filename.
    """

    sav_path: Path
    observation_date: str
    filters: list[str]
    intensities: list[float]
    label: str
    mc_runs: int = 0

    @property
    def intensities_array(self) -> np.ndarray:
        return np.array(self.intensities, dtype=float)

    @property
    def is_mc(self) -> bool:
        return self.mc_runs > 0


@dataclass(frozen=True)
class IDLResult:
    """Loaded IDL DEM solution."""

    logT: np.ndarray  # noqa: N815
    dem: np.ndarray  # (nT,) cm^-5 K^-1


@dataclass(frozen=True)
class IDLMCResult:
    """
    Loaded IDL DEM solution including all Monte Carlo runs.

    Attributes
    ----------
    logT : ndarray (nT,)
        log10 temperature grid.
    dem_base : ndarray (nT,)
        Base DEM — standalone base solution (loaded from base SAV).
    dem_mc : ndarray (n_runs, nT)
        Monte Carlo DEM realizations — columns 1..n_runs of dem_out.
    n_runs : int
        Number of MC realizations (excludes base).
    chisq_base : float
        Chi-square of the base DEM solution (run 0).
    chisq_mc : ndarray (n_runs,) or None
        Chi-square values for each MC realization.
        None if chi-square was not stored in the SAV file.
    """

    logT: np.ndarray  # noqa: N815
    dem_base: np.ndarray
    dem_mc: np.ndarray
    n_runs: int
    chisq_base: float = 0.0
    chisq_mc: np.ndarray | None = None


# ---------------------------------------------------------------------------
# Filename parser
# ---------------------------------------------------------------------------


def _parse_filter_token(token: str) -> tuple[str, float] | tuple[None, None]:
    for key in _FILTER_KEYS_SORTED:
        if token.startswith(key):
            remainder = token[len(key) :]
            try:
                return _FILTER_MAP[key], float(remainder)
            except ValueError:
                continue
    return None, None


def _parse_mc_token(token: str) -> int:
    """Parse 'MC1000' -> 1000. Returns 0 if not an MC token."""
    m = re.match(r"^MC(\d+)$", token, re.IGNORECASE)
    return int(m.group(1)) if m else 0


def parse_sav_filename(sav_path: str | Path) -> SavCase:
    """
    Parse an IDL DEM .sav filename and return a SavCase.

    Handles both base and MC filenames:
        xrt_IDL_dem_<DATE>_<Filter><I>_..._.sav
        xrt_IDL_dem_<DATE>_<Filter><I>_..._MC1000_.sav
    """
    sav_path = Path(sav_path)
    stem = sav_path.stem
    base = stem.removeprefix("xrt_IDL_dem_").strip("_")
    parts = [p for p in base.split("_") if p]

    if not parts:
        raise ValueError(f"Cannot parse filename: {sav_path.name}")

    # Date
    date_token = parts[0]
    m = re.match(r"(\d{4})(\d{2})(\d{2})T(\d{2})(\d{2})", date_token)

    if m:
        obs_date = f"{m.group(1)}-{m.group(2)}-{m.group(3)}T{m.group(4)}:{m.group(5)}"
    else:
        obs_date = date_token  # fall back to raw string

    filters: list[str] = []
    intensities: list[float] = []
    mc_runs = 0

    for token in parts[1:]:
        mc = _parse_mc_token(token)
        if mc > 0:
            mc_runs = mc
            continue
        fname, val = _parse_filter_token(token)
        if fname is not None:
            filters.append(fname)
            intensities.append(val)

    if not filters:
        raise ValueError(
            f"No filter/intensity pairs found in filename: {sav_path.name}"
        )

    return SavCase(
        sav_path=sav_path,
        observation_date=obs_date,
        filters=filters,
        intensities=intensities,
        label=date_token,
        mc_runs=mc_runs,
    )


# SAV load
def load_idl_sav(sav_path: str | Path) -> IDLResult:
    """Load base DEM only (column 0) from an IDL .sav file."""
    sav_path = Path(sav_path)
    data = readsav(str(sav_path), python_dict=True)

    logT = None
    for key in ("logt", "logT", "logT_out", "logt_out"):
        if key in data:
            logT = np.array(data[key]).ravel().astype(float)
            break
    if logT is None:
        raise KeyError(
            f"logT not found in {sav_path.name}. Keys: {sorted(data.keys())}"
        )

    nT = len(logT)
    dem = None
    for key in ("dem", "dem_out", "dem0"):
        if key in data:
            arr = np.array(data[key])
            if arr.ndim == 2:
                arr = arr[:, 0] if arr.shape[0] == nT else arr[0, :]
            dem = arr.ravel().astype(float)
            break
    if dem is None:
        raise KeyError(f"DEM not found in {sav_path.name}. Keys: {sorted(data.keys())}")

    return IDLResult(logT=logT, dem=dem)


def load_idl_mc_sav(sav_path: str | Path) -> IDLMCResult:
    """
    Load an IDL DEM .sav file containing Monte Carlo runs.

    Expected SAV structure:
        logT_out   : (nT,)
        dem_out    : (nT, nRuns+1)  or  (nRuns+1, nT)
        chisq_out  : (nRuns+1,)  or scalar  [optional]

    Returns IDLMCResult with:
        logT        (nT,)
        dem_base    (nT,)         — NOTE: this is dem_out[0], not the standalone base
        dem_mc      (n_runs, nT)
        n_runs      int
        chisq_base  float         — chi-square of run 0
        chisq_mc    (n_runs,) or None
    """
    sav_path = Path(sav_path)
    data = readsav(str(sav_path), python_dict=True)

    # logT
    logT = None
    for key in ("logt", "logT", "logT_out", "logt_out"):
        if key in data:
            logT = np.array(data[key]).ravel().astype(float)
            break
    if logT is None:
        raise KeyError(
            f"logT not found in {sav_path.name}. Keys: {sorted(data.keys())}"
        )

    nT = len(logT)

    # Load the complete DEM array.
    dem_all = None
    for key in ("dem_out", "dem", "dem0"):
        if key in data:
            dem_all = np.asarray(data[key], dtype=float).squeeze()
            break

    if dem_all is None:
        raise KeyError(
            f"DEM array not found in {sav_path.name}. "
            f"Keys: {sorted(data.keys())}"
        )

    if dem_all.ndim == 1:
        dem_all = dem_all[np.newaxis, :]
    elif dem_all.ndim != 2:
        raise ValueError(
            f"Expected a 1-D or 2-D DEM array, got shape {dem_all.shape}."
        )

    # Determine the number of MC runs recorded in the SAV file.
    saved_n_runs = 0
    if "monte_carlo_runs" in data:
        saved_n_runs = int(np.asarray(data["monte_carlo_runs"]).squeeze())

    # Normalize to (base + MC runs, temperature bins).
    expected_run_count = saved_n_runs + 1

    if dem_all.shape[0] == expected_run_count:
        runs_by_temperature = dem_all
    elif dem_all.shape[1] == expected_run_count:
        runs_by_temperature = dem_all.T
    else:
        raise ValueError(
            f"Cannot identify the run axis in DEM shape {dem_all.shape}; "
            f"expected one axis to have length {expected_run_count}."
        )

    n_temperature_bins = runs_by_temperature.shape[1]

    # IDL may save temperature-bin edges rather than bin centers.
    if logT.size == n_temperature_bins + 1:
        logT = 0.5 * (logT[:-1] + logT[1:])
    elif logT.size != n_temperature_bins:
        raise ValueError(
            f"Temperature-grid length {logT.size} is incompatible with "
            f"{n_temperature_bins} DEM temperature bins."
        )

    dem_base = runs_by_temperature[0].copy()
    dem_mc = runs_by_temperature[1:].copy()
    n_runs = dem_mc.shape[0]



    # chi-square — try to load, handle gracefully if absent or wrong shape
    chisq_base = 0.0
    chisq_mc = None

    for key in ("chisq_out", "chisq", "chi2_out", "chi2"):
        if key in data:
            raw = np.array(data[key]).ravel().astype(float)

            if raw.size == 1:
                # Scalar — only base chi-square was saved
                chisq_base = float(raw[0])
                chisq_mc = None

            elif raw.size == n_runs + 1:
                # Full array: index 0 = base, 1..n_runs = MC
                chisq_base = float(raw[0])
                chisq_mc = raw[1:].copy()

            elif raw.size == n_runs:
                # MC only (no base stored)
                chisq_base = 0.0
                chisq_mc = raw.copy()

            else:
                # Unexpected size — store what we can
                chisq_base = float(raw[0]) if raw.size > 0 else 0.0
                chisq_mc = raw[1:] if raw.size > 1 else None

            break  # found a chi-square key, stop looking

    return IDLMCResult(
        logT=logT,
        dem_base=dem_base,
        dem_mc=dem_mc,
        n_runs=n_runs,
        chisq_base=chisq_base,
        chisq_mc=chisq_mc,
    )


# ---------------------------------------------------------------------------
# Directory discovery
# ---------------------------------------------------------------------------


def discover_cases(data_dir: str | Path) -> list[SavCase]:
    """Find all xrt_IDL_dem_*.sav files in data_dir and parse each."""
    data_dir = Path(data_dir)
    cases: list[SavCase] = []
    for sav_file in sorted(data_dir.glob("xrt_IDL_dem_*.sav")):
        try:
            cases.append(parse_sav_filename(sav_file))
        except ValueError as exc:  # noqa: PERF203
            print(f"  Warning: skipping {sav_file.name} — {exc}")
    return cases


def discover_mc_cases(mc_dir: str | Path) -> list[SavCase]:
    """Find all MC .sav files (those with MC<N> token in filename)."""
    cases = discover_cases(mc_dir)
    return [c for c in cases if c.is_mc]

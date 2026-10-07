"""
Scientific Validation: IDL xrt_dem_iterative2.pro vs XRTpy XRTDEMIterative
============================================================================

These tests compare the XRTpy DEM solver output against reference solutions
produced by the IDL routine xrt_dem_iterative2.pro (SolarSoft).

How to add a new case
---------------------
Drop a new .sav file into ``data/validation/`` using this naming convention::

    xrt_IDL_dem_<YYYYMMDDTHHMI>_<Filter1><Intensity1>_<Filter2><Intensity2>_..._.sav

Example::

    xrt_IDL_dem_20071213T0401_Bemed603.875886_Bethin150.921435_Alpoly2412.34_.sav

The test suite discovers and runs all matching files automatically.
No code changes required.

Tolerances (Standard tier)
--------------------------
    mean |Δlog10(DEM)| < 0.20 dex
    max  |Δlog10(DEM)| < 0.50 dex
    peak logT difference < 0.15 dex (~1 bin)
"""

from pathlib import Path

import numpy as np
import pytest
from utils_sav_io import (
    IDLResult,
    SavCase,
    discover_cases,
    load_idl_sav,
)

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ---------------------------------------------------------------------------
# Tolerances
# ---------------------------------------------------------------------------
MEAN_DEX_TOL = 0.20
MAX_DEX_TOL = 0.50
PEAK_LOGT_TOL = 0.15
DEM_FLOOR = 1e10  # cm^-5 K^-1 — bins below this are ignored


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _log10_safe(arr: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(arr, floor))


def _valid_mask(dem_idl: np.ndarray, dem_xrt: np.ndarray) -> np.ndarray:
    """Bins where at least one DEM is physically meaningful."""
    return (dem_idl > DEM_FLOOR) | (dem_xrt > DEM_FLOOR)


# ---------------------------------------------------------------------------
# Case discovery
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data" / "validation" / "base"


def _collect_cases() -> list[SavCase]:
    if not DATA_DIR.exists():
        return []
    return discover_cases(DATA_DIR)


CASES = _collect_cases()


def _case_id(case: SavCase) -> str:
    return case.label


# ---------------------------------------------------------------------------
# Session-scoped fixture: solve XRTpy once per case, reused across all tests
# ---------------------------------------------------------------------------


@pytest.fixture(scope="session", params=CASES, ids=_case_id)
def solved(request) -> tuple[SavCase, IDLResult, XRTDEMIterative]:
    """
    Returns (case, idl_result, xrtpy_solver) for one .sav case.
    Skips if the .sav file is missing.
    XRTpy is solved once and shared across all tests for that case.
    """
    case: SavCase = request.param

    if not case.sav_path.exists():
        pytest.skip(f"SAV file not found: {case.sav_path}")

    idl = load_idl_sav(case.sav_path)

    responses = generate_temperature_responses(
        case.filters,
        idl.observation_date,
    )

    solver = XRTDEMIterative(
        observed_channel=case.filters,
        observed_intensities=idl.observed_intensities,
        temperature_responses=responses,
        monte_carlo_runs=0,
    )
    solver.solve()

    return case, idl, solver


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


def test_logt_grids_are_consistent(solved):
    """IDL and XRTpy must use the same logT grid."""
    case, idl, xrtpy = solved
    assert idl.logT.shape == xrtpy.logT.shape, (
        f"[{case.label}] Grid size mismatch: IDL={idl.logT.shape}, XRTpy={xrtpy.logT.shape}"
    )
    np.testing.assert_allclose(
        idl.logT,
        xrtpy.logT,
        atol=1e-6,
        err_msg=f"[{case.label}] logT grids differ",
    )


def test_mean_log10_dem_difference(solved):
    """Mean |Δlog10(DEM)| must be < MEAN_DEX_TOL across valid bins."""
    case, idl, xrtpy = solved
    mask = _valid_mask(idl.dem, xrtpy.dem)
    assert mask.sum() >= 5, f"[{case.label}] Too few valid bins"

    diff = np.abs(_log10_safe(xrtpy.dem[mask]) - _log10_safe(idl.dem[mask]))
    mean_diff = float(np.mean(diff))

    print(
        f"\n  [{case.label}]  Mean |Δlog10(DEM)| = {mean_diff:.4f} dex  (tol={MEAN_DEX_TOL})"
    )
    assert mean_diff < MEAN_DEX_TOL, (
        f"[{case.label}] Mean Δ = {mean_diff:.3f} dex > {MEAN_DEX_TOL}"
    )


def test_max_log10_dem_difference(solved):
    """Max |Δlog10(DEM)| must be < MAX_DEX_TOL across valid bins."""
    case, idl, xrtpy = solved
    mask = _valid_mask(idl.dem, xrtpy.dem)

    diff = np.abs(_log10_safe(xrtpy.dem[mask]) - _log10_safe(idl.dem[mask]))
    max_diff = float(np.max(diff))

    print(
        f"\n  [{case.label}]  Max |Δlog10(DEM)| = {max_diff:.4f} dex  (tol={MAX_DEX_TOL})"
    )
    assert max_diff < MAX_DEX_TOL, (
        f"[{case.label}] Max Δ = {max_diff:.3f} dex > {MAX_DEX_TOL}"
    )


def test_peak_temperature_agreement(solved):
    """Peak logT must agree within PEAK_LOGT_TOL."""
    case, idl, xrtpy = solved
    pk_idl = idl.logT[np.argmax(idl.dem)]
    pk_xrt = xrtpy.logT[np.argmax(xrtpy.dem)]
    diff = abs(pk_xrt - pk_idl)

    print(f"\n  [{case.label}]  IDL={pk_idl:.2f}  XRTpy={pk_xrt:.2f}  Δ={diff:.3f}")
    assert diff < PEAK_LOGT_TOL, (
        f"[{case.label}] Peak logT Δ={diff:.3f} > {PEAK_LOGT_TOL} "
        f"(IDL={pk_idl:.2f}, XRTpy={pk_xrt:.2f})"
    )


def test_modeled_intensities_are_finite_and_positive(solved):
    """Modeled intensities must be finite and non-negative."""
    case, _, xrtpy = solved
    assert np.all(np.isfinite(xrtpy.modeled_intensities)), (
        f"[{case.label}] Non-finite modeled intensity"
    )
    assert np.all(xrtpy.modeled_intensities >= 0.0), (
        f"[{case.label}] Negative modeled intensity"
    )


def test_modeled_intensities_order_of_magnitude(solved):
    """Modeled intensities must be within 2 dex of observed."""
    case, _, xrtpy = solved
    ratio = xrtpy.modeled_intensities / case.intensities_array
    log_ratio = np.log10(np.maximum(ratio, 1e-99))

    print(f"\n  [{case.label}]  log10(mod/obs):")
    for f, lr in zip(case.filters, log_ratio):
        print(f"    {f:<22} {lr:+.3f}{'  ← !' if abs(lr) > 1.0 else ''}")

    assert np.all(np.abs(log_ratio) < 2.0), (
        f"[{case.label}] Modeled intensity >2 dex from observed.\n"
        f"  Filters:    {case.filters}\n"
        f"  log10(M/O): {log_ratio.round(3)}"
    )


def test_chisq_is_finite(solved):
    """Chi-square must be finite."""
    case, _, xrtpy = solved
    assert np.isfinite(xrtpy.chisq), f"[{case.label}] χ² is not finite: {xrtpy.chisq}"


def test_chisq_is_reasonable(solved):
    """Reduced chi-square must be < 10."""
    case, _, xrtpy = solved
    n_dof = max(1, len(case.filters) - xrtpy.n_spl)
    reduced = xrtpy.chisq / n_dof
    print(
        f"\n  [{case.label}]  χ²={xrtpy.chisq:.2f}  reduced χ²={reduced:.2f}  dof={n_dof}"
    )
    assert reduced < 10.0, f"[{case.label}] Reduced χ² = {reduced:.2f} > 10"


def test_diagnostic_print_full_comparison(solved):
    """
    Always passes.  Prints the full per-bin table and per-filter breakdown.
    Use ``pytest -s`` to see output.
    """
    case, idl, xrtpy = solved
    log_idl = _log10_safe(idl.dem)
    log_xrt = _log10_safe(xrtpy.dem)
    mask = _valid_mask(idl.dem, xrtpy.dem)
    delta = log_xrt - log_idl

    mean_d = float(np.mean(np.abs(delta[mask])))
    max_d = float(np.max(np.abs(delta[mask])))

    sep = "=" * 65
    print(f"\n{sep}")
    print(f"  IDL vs XRTpy  |  Case: {case.label}")
    print(f"  Date:     {case.observation_date}")
    print(f"  Filters:  {case.filters}")
    print(f"{sep}")
    print(f"  {'logT':>6}  {'IDL':>10}  {'XRTpy':>10}  {'Δ(dex)':>9}  valid")
    print(f"  {'-' * 55}")
    for i, lt in enumerate(idl.logT):
        d = delta[i] if mask[i] else float("nan")
        flag = " *" if mask[i] and abs(d) > MAX_DEX_TOL else "  "
        print(
            f"  {lt:>6.2f}  {log_idl[i]:>10.3f}  {log_xrt[i]:>10.3f}  "
            f"{d:>+9.3f}  {'yes' if mask[i] else 'no '}{flag}"
        )
    print(f"  {'-' * 55}")
    print(f"  Mean |Δ| = {mean_d:.4f} dex     Max |Δ| = {max_d:.4f} dex")
    print(f"  XRTpy χ² = {xrtpy.chisq:.2f}")
    print(f"  IDL   peak logT = {idl.logT[np.argmax(idl.dem)]:.2f}")
    print(f"  XRTpy peak logT = {xrtpy.logT[np.argmax(xrtpy.dem)]:.2f}")
    print()
    print(f"  {'Filter':<22} {'Observed':>10} {'Modeled':>10} {'log10(M/O)':>11}")
    print(f"  {'-' * 55}")
    for f, obs, mod in zip(
        case.filters, case.intensities_array, xrtpy.modeled_intensities
    ):
        lr = np.log10(max(mod / obs, 1e-99))
        print(
            f"  {f:<22} {obs:>10.3f} {mod:>10.3f} {lr:>+11.3f}{'  ←!' if abs(lr) > 1.0 else ''}"
        )
    print(f"{sep}\n")


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

from xrtpy.xrt_dem_iterative.utils_sav_io import IDLResult, SavCase, discover_cases, load_idl_sav

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

DATA_DIR = Path(__file__).parent / "data" / "validation"


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

    responses = generate_temperature_responses(case.filters, case.observation_date)
    solver = XRTDEMIterative(
        observed_channel=case.filters,
        observed_intensities=case.intensities_array,
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

    print(f"\n  [{case.label}]  Mean |Δlog10(DEM)| = {mean_diff:.4f} dex  (tol={MEAN_DEX_TOL})")
    assert mean_diff < MEAN_DEX_TOL, (
        f"[{case.label}] Mean Δ = {mean_diff:.3f} dex > {MEAN_DEX_TOL}"
    )


def test_max_log10_dem_difference(solved):
    """Max |Δlog10(DEM)| must be < MAX_DEX_TOL across valid bins."""
    case, idl, xrtpy = solved
    mask = _valid_mask(idl.dem, xrtpy.dem)

    diff = np.abs(_log10_safe(xrtpy.dem[mask]) - _log10_safe(idl.dem[mask]))
    max_diff = float(np.max(diff))

    print(f"\n  [{case.label}]  Max |Δlog10(DEM)| = {max_diff:.4f} dex  (tol={MAX_DEX_TOL})")
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
    assert np.isfinite(xrtpy.chisq), (
        f"[{case.label}] χ² is not finite: {xrtpy.chisq}"
    )


def test_chisq_is_reasonable(solved):
    """Reduced chi-square must be < 10."""
    case, _, xrtpy = solved
    n_dof = max(1, len(case.filters) - xrtpy.n_spl)
    reduced = xrtpy.chisq / n_dof
    print(f"\n  [{case.label}]  χ²={xrtpy.chisq:.2f}  reduced χ²={reduced:.2f}  dof={n_dof}")
    assert reduced < 10.0, (
        f"[{case.label}] Reduced χ² = {reduced:.2f} > 10"
    )


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
    print(f"  {'-'*55}")
    for i, lt in enumerate(idl.logT):
        d = delta[i] if mask[i] else float("nan")
        flag = " *" if mask[i] and abs(d) > MAX_DEX_TOL else "  "
        print(
            f"  {lt:>6.2f}  {log_idl[i]:>10.3f}  {log_xrt[i]:>10.3f}  "
            f"{d:>+9.3f}  {'yes' if mask[i] else 'no '}{flag}"
        )
    print(f"  {'-'*55}")
    print(f"  Mean |Δ| = {mean_d:.4f} dex     Max |Δ| = {max_d:.4f} dex")
    print(f"  XRTpy χ² = {xrtpy.chisq:.2f}")
    print(f"  IDL   peak logT = {idl.logT[np.argmax(idl.dem)]:.2f}")
    print(f"  XRTpy peak logT = {xrtpy.logT[np.argmax(xrtpy.dem)]:.2f}")
    print()
    print(f"  {'Filter':<22} {'Observed':>10} {'Modeled':>10} {'log10(M/O)':>11}")
    print(f"  {'-'*55}")
    for f, obs, mod in zip(case.filters, case.intensities_array, xrtpy.modeled_intensities):
        lr = np.log10(max(mod / obs, 1e-99))
        print(f"  {f:<22} {obs:>10.3f} {mod:>10.3f} {lr:>+11.3f}{'  ←!' if abs(lr) > 1.0 else ''}")
    print(f"{sep}\n")


# """
# Scientific Validation: IDL xrt_dem_iterative2.pro vs XRTpy XRTDEMIterative
# ============================================================================

# These tests compare the XRTpy DEM solver output against reference solutions
# produced by the IDL routine xrt_dem_iterative2.pro (SolarSoft).

# Because DEM inversion is ill-posed and the two implementations use different
# optimizers (IDL MPFIT vs Python lmfit) and spline libraries, exact bin-by-bin
# agreement is not expected. Instead, we verify:

#     1. log10(DEM) mean absolute difference < tolerance (in dex)
#     2. log10(DEM) max absolute difference  < tolerance (in dex)
#     3. Peak temperature agreement within one logT bin (0.1 dex)
#     4. Modeled intensities agree with observed within chi-square tolerance

# Tolerances used (Standard tier):
#     mean |log10(DEM)| < 0.20 dex
#     max  |log10(DEM)| < 0.50 dex
#     peak logT difference < 0.15 (within ~1 bin)
# """

# from pathlib import Path

# import numpy as np
# import pytest
# from scipy.io import readsav

# from xrtpy.response.tools import generate_temperature_responses
# from xrtpy.xrt_dem_iterative import XRTDEMIterative


# DATA_DIR = Path(__file__).parent / "data" / "validation"

# # Tolerance tier — Standard
# MEAN_DEX_TOL = 0.20   # mean |log10(DEM)| across valid bins
# MAX_DEX_TOL  = 0.50   # max  |log10(DEM)| across valid bins
# PEAK_LOGT_TOL = 0.15  # peak temperature agreement [log10 K]

# # Mask threshold — ignore bins where both DEMs are essentially zero
# DEM_FLOOR = 1e10  # [cm^-5 K^-1]  below this both are noise


# def _log10_safe(dem, floor=1e-99):
#     """log10 with a floor to avoid -inf."""
#     return np.log10(np.maximum(dem, floor))


# def _valid_mask(dem_idl, dem_xrtpy, floor=DEM_FLOOR):
#     """
#     Boolean mask of bins where at least one DEM is above the floor.
#     Excludes essentially-zero tails from the comparison.
#     """
#     return (dem_idl > floor) | (dem_xrtpy > floor)


# def _load_idl_sav(sav_path):
#     """
#     Load IDL .sav file and return (logT, dem_base) as 1D float arrays.
#     """
#     data = readsav(str(sav_path), python_dict=True)

#     # logT — try common key names
#     for key in ("logt", "logT", "logT_out", "logt_out"):
#         if key in data:
#             logT = np.array(data[key]).ravel().astype(float)
#             break
#     else:
#         raise KeyError(f"logT key not found. Available keys: {list(data.keys())}")

#     # dem — base DEM (1D or first column of 2D)
#     for key in ("dem", "dem_out", "dem0"):
#         if key in data:
#             dem = np.array(data[key])
#             if dem.ndim == 2:
#                 dem = dem[:, 0]   # first column = base
#             dem = dem.ravel().astype(float)
#             break
#     else:
#         raise KeyError(f"DEM key not found. Available keys: {list(data.keys())}")

#     return logT, dem


# # ---------------------------------------------------------------------------
# # Case 1: 2008-01-04  (5-filter, no MC)
# # ---------------------------------------------------------------------------

# class TestIDLvsXRTpy_20080104:


#     SAV_FILE = DATA_DIR / "xrt_dem_output_20071213T0401_NOMC.sav"#"xrt_IDL_dem_output_20080104.sav"

    FILTERS = ["Be-med", "Al-mesh", "Ti-poly", "Al-poly", "Be-thin"]
    INTENSITIES = np.array([
        234.283365,  # Be-med
        183.711876,  # Al-mesh
        45.931438,  # Ti-poly
        91.745329,  # Al-poly
        5.755926,  # Be-thin
        ], dtype=float)
    OBSERVATION_DATE = "2008-01-04T11:04:26"
    
#     FILTERS= ["be-med","Be-thin","Al-poly", "Al-poly/Ti-poly","Ti-poly","Al-thick"]
#     INTENSITIES = [603.875886,150.921435,2412.340960, 301.354389 ,603.100596,2.519851]
#     OBSERVATION_DATE = "2007-12-13T04:01"

#     @pytest.fixture(scope="class")
#     def idl(self):
#         """Load IDL reference solution."""
#         if not self.SAV_FILE.exists():
#             pytest.skip(f"IDL .sav file not found: {self.SAV_FILE}")
#         logT, dem = _load_idl_sav(self.SAV_FILE)
#         return {"logT": logT, "dem": dem}

#     @pytest.fixture(scope="class")
#     def xrtpy(self):
#         """Run XRTpy solver with the same inputs."""
#         responses = generate_temperature_responses(
#             self.FILTERS, self.OBSERVATION_DATE
#         )
#         solver = XRTDEMIterative(
#             observed_channel=self.FILTERS,
#             observed_intensities=self.INTENSITIES,
#             temperature_responses=responses,
#         )
#         solver.solve()
#         return solver

#     def test_logt_grids_are_consistent(self, idl, xrtpy):
#         """IDL and XRTpy must use the same logT grid."""
#         assert idl["logT"].shape == xrtpy.logT.shape, (
#             f"Grid size mismatch: IDL={idl['logT'].shape}, XRTpy={xrtpy.logT.shape}"
#         )
#         np.testing.assert_allclose(
#             idl["logT"], xrtpy.logT, atol=1e-6,
#             err_msg="logT grids differ between IDL and XRTpy"
#         )

#     def test_mean_log10_dem_difference(self, idl, xrtpy):
#         """
#         Mean |log10(DEM)| across valid bins must be < MEAN_DEX_TOL.
#         Valid bins are those where at least one DEM exceeds DEM_FLOOR.
#         """
#         mask = _valid_mask(idl["dem"], xrtpy.dem)
#         assert mask.sum() >= 5, "Too few valid bins for comparison"

#         log_idl  = _log10_safe(idl["dem"][mask])
#         log_xrt  = _log10_safe(xrtpy.dem[mask])
#         mean_diff = np.mean(np.abs(log_xrt - log_idl))

#         print(f"\n  Mean |Δlog10(DEM)| = {mean_diff:.4f} dex  (tolerance={MEAN_DEX_TOL})")
#         assert mean_diff < MEAN_DEX_TOL, (
#             f"Mean log10(DEM) difference {mean_diff:.3f} dex exceeds "
#             f"tolerance {MEAN_DEX_TOL} dex"
#         )

#     def test_max_log10_dem_difference(self, idl, xrtpy):
#         """
#         Max |Δlog10(DEM)| across valid bins must be < MAX_DEX_TOL.
#         """
#         mask = _valid_mask(idl["dem"], xrtpy.dem)

#         log_idl  = _log10_safe(idl["dem"][mask])
#         log_xrt  = _log10_safe(xrtpy.dem[mask])
#         max_diff  = np.max(np.abs(log_xrt - log_idl))

#         print(f"\n  Max  |Δlog10(DEM)| = {max_diff:.4f} dex  (tolerance={MAX_DEX_TOL})")
#         assert max_diff < MAX_DEX_TOL, (
#             f"Max log10(DEM) difference {max_diff:.3f} dex exceeds "
#             f"tolerance {MAX_DEX_TOL} dex"
#         )

#     def test_peak_temperature_agreement(self, idl, xrtpy):
#         """
#         The temperature at which DEM peaks must agree within PEAK_LOGT_TOL.
#         """
#         peak_logT_idl  = idl["logT"][np.argmax(idl["dem"])]
#         peak_logT_xrt  = xrtpy.logT[np.argmax(xrtpy.dem)]
#         diff = abs(peak_logT_xrt - peak_logT_idl)

#         print(f"\n  IDL  peak logT = {peak_logT_idl:.2f}")
#         print(f"  XRTpy peak logT = {peak_logT_xrt:.2f}")
#         print(f"  Difference      = {diff:.3f}  (tolerance={PEAK_LOGT_TOL})")

#         assert diff < PEAK_LOGT_TOL, (
#             f"Peak logT differs by {diff:.3f} dex "
#             f"(IDL={peak_logT_idl:.2f}, XRTpy={peak_logT_xrt:.2f})"
#         )


#     def test_modeled_intensities_are_finite_and_positive(self, xrtpy):
#         """Modeled intensities must be finite and non-negative."""
#         assert np.all(np.isfinite(xrtpy.modeled_intensities))
#         assert np.all(xrtpy.modeled_intensities >= 0.0)

#     def test_modeled_intensities_order_of_magnitude(self, xrtpy):
#         """
#         Modeled intensities should be within 2 orders of magnitude
#         of observed intensities (basic physical sanity check).
#         """
#         ratio = xrtpy.modeled_intensities / self.INTENSITIES
#         log_ratio = np.log10(np.maximum(ratio, 1e-99))

#         print(f"\n  log10(modeled/observed) per filter: {log_ratio.round(3)}")
#         assert np.all(np.abs(log_ratio) < 2.0), (
#             f"Modeled intensities differ from observed by more than "
#             f"2 orders of magnitude: log10(ratio)={log_ratio}"
#         )

#     def test_chisq_is_finite(self, xrtpy):
#         """Chi-square of the base fit must be finite."""
#         assert np.isfinite(xrtpy.chisq), f"Chi-square is not finite: {xrtpy.chisq}"

#     def test_chisq_is_reasonable(self, xrtpy):
#         """
#         Reduced chi-square should be < 10 for a reasonable fit.
#         (n_filters=5, n_knots=4 → dof=1, so chi-square itself should be small.)
#         """
#         n_filters = len(self.FILTERS)
#         reduced_chisq = xrtpy.chisq / max(1, n_filters - xrtpy.n_spl)
#         print(f"\n  chi-square         = {xrtpy.chisq:.4f}")
#         print(f"  reduced chi-square = {reduced_chisq:.4f}")
#         assert reduced_chisq < 10.0, (
#             f"Reduced chi-square {reduced_chisq:.2f} is unexpectedly large"
#         )
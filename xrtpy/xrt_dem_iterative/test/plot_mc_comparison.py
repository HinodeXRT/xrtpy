# """
# plot_mc_comparison.py
# =====================
# Four-panel Monte Carlo comparison: IDL vs XRTpy DEM ensembles.

# Usage
# -----
#     python plot_mc_comparison.py --case 20071213T0401
#     python plot_mc_comparison.py --all
#     python plot_mc_comparison.py --all --save

# The script expects:
#     data/validation/base/        xrt_IDL_dem_<LABEL>_..._.sav
#     data/validation/monte_carlo/ xrt_IDL_dem_<LABEL>_..._MC<N>_.sav

# Plots produced per case:
#     plots/<LABEL>/mc_plot1_idlmc_xrtpy_base.png
#     plots/<LABEL>/mc_plot2_idl_base_xrtpy_mc.png
#     plots/<LABEL>/mc_plot3_both_ensembles.png
#     plots/<LABEL>/mc_plot4_combined.png        ← shareable 3-panel figure
# """

# import argparse
# import sys
# import warnings
# from pathlib import Path

# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# import numpy as np

# from utils_sav_io import load_idl_sav

# from utils_sav_io import (
#     SavCase,
#     IDLMCResult,
#     discover_mc_cases,
#     load_idl_mc_sav,
#     parse_sav_filename,
# )

# from xrtpy.response.tools import generate_temperature_responses
# from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # ---------------------------------------------------------------------------
# # Config
# # ---------------------------------------------------------------------------
# BASE_DIR = Path(__file__).parent / "data" / "validation" / "base"
# MC_DIR   = Path(__file__).parent / "data" / "validation" / "monte_carlo"
# PLOT_DIR = Path(__file__).parent / "plots"

# # Visual settings
# IDL_COLOR    = "#BA7517"   # amber — IDL
# XRTPY_COLOR  = "#185FA5"   # blue  — XRTpy
# IDL_ALPHA    = 0.12        # transparency for ensemble curves
# XRTPY_ALPHA  = 0.12
# BASE_LW      = 2.5         # linewidth for the single base DEM
# ENS_LW       = 0.6         # linewidth for ensemble curves
# DEM_FLOOR    = 1e-99


# def _log10(dem: np.ndarray) -> np.ndarray:
#     return np.log10(np.maximum(dem, DEM_FLOOR))


# def _axis_style(ax, title: str, xlim=(5.4, 8.1)):
#     ax.set_xlabel("log₁₀ T  [K]", fontsize=10)
#     ax.set_ylabel("log₁₀ DEM  [cm⁻⁵ K⁻¹]", fontsize=10)
#     ax.set_title(title, fontsize=11)
#     ax.set_xlim(*xlim)
#     ax.grid(True, alpha=0.22)


# # ---------------------------------------------------------------------------
# # Solver
# # ---------------------------------------------------------------------------

# def run_xrtpy(case: SavCase, n_mc: int) -> XRTDEMIterative:
#     """Run XRTpy solver with the same inputs as the IDL case."""
#     with warnings.catch_warnings():
#         warnings.simplefilter("ignore")
#         responses = generate_temperature_responses(case.filters, case.observation_date)
#         solver = XRTDEMIterative(
#             observed_channel=case.filters,
#             observed_intensities=case.intensities_array,
#             temperature_responses=responses,
#             monte_carlo_runs=n_mc,
#         )
#         solver.solve()
#     return solver


# # ---------------------------------------------------------------------------
# # The 3 plot functions + combined
# # ---------------------------------------------------------------------------

# def plot1_idlmc_xrtpy_base(
#     idl: IDLMCResult,
#     xrtpy_base: np.ndarray,
#     logT: np.ndarray,
#     case: SavCase,
#     ax=None,
# ) -> plt.Figure | None:
#     """
#     Plot 1: IDL 1000 MC ensemble (thin amber) + XRTpy base DEM (thick blue).
#     """
#     standalone = ax is None
#     if standalone:
#         fig, ax = plt.subplots(figsize=(9, 5.5))

#     # IDL MC ensemble
#     for i in range(idl.n_runs):
#         ax.step(idl.logT, _log10(idl.dem_mc[i]), where="mid",
#                 color=IDL_COLOR, alpha=IDL_ALPHA, linewidth=ENS_LW)

#     # IDL base
#     ax.step(idl.logT, _log10(idl.dem_base), where="mid",
#             color=IDL_COLOR, linewidth=BASE_LW, label=f"IDL base DEM")

#     # XRTpy base (thick blue dashed)
#     ax.step(logT, _log10(xrtpy_base), where="mid",
#             color=XRTPY_COLOR, linewidth=BASE_LW, linestyle="--",
#             label="XRTpy base DEM")

#     _axis_style(ax, f"Plot 1 — IDL {idl.n_runs} MC + XRTpy base")
#     ax.legend(fontsize=9)

#     if standalone:
#         fig.suptitle(f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
#                      fontsize=10, y=1.01)
#         fig.tight_layout()
#         return fig
#     return None


# def plot2_idl_base_xrtpy_mc(
#     idl_base: np.ndarray,
#     idl_logT: np.ndarray,
#     xrtpy: XRTDEMIterative,
#     case: SavCase,
#     ax=None,
# ) -> plt.Figure | None:
#     """
#     Plot 2: IDL base DEM (thick amber) + XRTpy 1000 MC ensemble (thin blue).
#     """
#     standalone = ax is None
#     if standalone:
#         fig, ax = plt.subplots(figsize=(9, 5.5))

#     n_mc = xrtpy.mc_dem.shape[0] - 1   # excludes row 0 (base)

#     # XRTpy MC ensemble
#     for i in range(1, n_mc + 1):
#         ax.step(xrtpy.logT, _log10(xrtpy.mc_dem[i]), where="mid",
#                 color=XRTPY_COLOR, alpha=XRTPY_ALPHA, linewidth=ENS_LW)

#     # XRTpy base
#     ax.step(xrtpy.logT, _log10(xrtpy.dem), where="mid",
#             color=XRTPY_COLOR, linewidth=BASE_LW, linestyle="--",
#             label="XRTpy base DEM")

#     # IDL base (thick amber solid)
#     ax.step(idl_logT, _log10(idl_base), where="mid",
#             color=IDL_COLOR, linewidth=BASE_LW,
#             label="IDL base DEM")

#     _axis_style(ax, f"Plot 2 — IDL base + XRTpy {n_mc} MC")
#     ax.legend(fontsize=9)

#     if standalone:
#         fig.suptitle(f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
#                      fontsize=10, y=1.01)
#         fig.tight_layout()
#         return fig
#     return None


# def plot3_both_ensembles(
#     idl: IDLMCResult,
#     xrtpy: XRTDEMIterative,
#     case: SavCase,
#     idl_true_base=None,
#     ax=None,
# ) -> plt.Figure | None:
#     if idl_true_base is None:
#         idl_true_base = idl.dem_base
#     """
#     Plot 3: IDL 1000 MC (amber) + XRTpy 1000 MC (blue) both overlaid.
#     """
#     standalone = ax is None
#     if standalone:
#         fig, ax = plt.subplots(figsize=(9, 5.5))

#     n_xrt = xrtpy.mc_dem.shape[0] - 1

#     # IDL ensemble
#     for i in range(idl.n_runs):
#         ax.step(idl.logT, _log10(idl.dem_mc[i]), where="mid",
#                 color=IDL_COLOR, alpha=IDL_ALPHA, linewidth=ENS_LW)

#     # XRTpy ensemble
#     for i in range(1, n_xrt + 1):
#         ax.step(xrtpy.logT, _log10(xrtpy.mc_dem[i]), where="mid",
#                 color=XRTPY_COLOR, alpha=XRTPY_ALPHA, linewidth=ENS_LW)

#     # Base DEMs (thick, on top)
#     ax.step(idl.logT, _log10(idl_true_base), where="mid",
#             color=IDL_COLOR, linewidth=BASE_LW, label=f"IDL base  ({idl.n_runs} MC shown)")
#     ax.step(xrtpy.logT, _log10(xrtpy.dem), where="mid",
#             color=XRTPY_COLOR, linewidth=BASE_LW, linestyle="--",
#             label=f"XRTpy base ({n_xrt} MC shown)")

#     _axis_style(ax, f"Plot 3 — Both ensembles overlaid\n"
#                     f"IDL {idl.n_runs} MC (amber) vs XRTpy {n_xrt} MC (blue)")
#     ax.legend(fontsize=9)

#     if standalone:
#         fig.suptitle(f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
#                      fontsize=10, y=1.01)
#         fig.tight_layout()
#         return fig
#     return None

# def plot4_combined(
#     idl: IDLMCResult,
#     xrtpy: XRTDEMIterative,
#     case: SavCase,
#     idl_true_base=None,
# ) -> plt.Figure:
#     if idl_true_base is None:
#         idl_true_base = idl.dem_base
#     """
#     Plot 4: All three panels stacked vertically — the shareable team figure.
#     """
#     fig = plt.figure(figsize=(10, 17))

#     fig.suptitle(
#         f"IDL vs XRTpy — Monte Carlo DEM Comparison\n"
#         f"Case: {case.observation_date}   "
#         f"Filters: {', '.join(case.filters)}\n"
#         f"IDL: {idl.n_runs} MC runs   |   "
#         f"XRTpy: {xrtpy.mc_dem.shape[0]-1} MC runs",
#         fontsize=12, fontweight="bold", y=0.99,
#     )

#     gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.42)

#     ax1 = fig.add_subplot(gs[0])
#     ax2 = fig.add_subplot(gs[1])
#     ax3 = fig.add_subplot(gs[2])

#     # plot1_idlmc_xrtpy_base(idl, xrtpy.dem, xrtpy.logT, case, ax=ax1)
#     # plot2_idl_base_xrtpy_mc(idl.dem_base, idl.logT, xrtpy, case, ax=ax2)
#     # plot3_both_ensembles(idl, xrtpy, case, ax=ax3)
#     plot1_idlmc_xrtpy_base(idl, xrtpy.dem, xrtpy.logT, case, ax=ax1)
#     plot2_idl_base_xrtpy_mc(idl_true_base, idl.logT, xrtpy, case, ax=ax2)
#     plot3_both_ensembles(idl, xrtpy, case, idl_true_base=idl_true_base, ax=ax3)

#     return fig


# # ---------------------------------------------------------------------------
# # Main
# # ---------------------------------------------------------------------------
# def process_case(mc_case: SavCase, save: bool) -> None:
#     """Load data, run XRTpy, generate all 4 plots for one case."""
#     print(f"\n{'='*60}")
#     print(f"  Case: {mc_case.label}  ({mc_case.observation_date})")
#     print(f"  Filters: {mc_case.filters}")
#     print(f"  IDL MC runs: {mc_case.mc_runs}")

#     # Load IDL MC results
#     print(f"  Loading MC SAV: {mc_case.sav_path.name}")
#     idl = load_idl_mc_sav(mc_case.sav_path)
#     print(f"  IDL dem_mc shape: {idl.dem_mc.shape}  (n_runs x nT)")

#     # Load the TRUE base DEM from the standalone base SAV
#     # (dem_out[0] from MC run is NOT the same as the standalone base solution)
#     base_sav_dir = MC_DIR.parent / "base"
#     base_matches = [
#         f for f in base_sav_dir.glob(f"xrt_IDL_dem_{mc_case.label}_*.sav")
#         if "MC" not in f.name
#     ]

#     if base_matches:
#         idl_base_result = load_idl_sav(base_matches[0])
#         idl_true_base = idl_base_result.dem
#         print(f"  Loaded true base DEM from: {base_matches[0].name}")
#     else:
#         print(f"  WARNING: no base SAV found for {mc_case.label} — falling back to MC row 0")
#         idl_true_base = idl.dem_base

#     # Run XRTpy with same number of MC runs
#     print(f"  Running XRTpy with {mc_case.mc_runs} MC runs...")
#     xrtpy = run_xrtpy(mc_case, n_mc=mc_case.mc_runs)
#     print(f"  XRTpy mc_dem shape: {xrtpy.mc_dem.shape}")
#     print(f"  XRTpy base χ²: {xrtpy.chisq:.2f}")

#     # Output folder
#     case_plot_dir = PLOT_DIR / mc_case.label
#     case_plot_dir.mkdir(parents=True, exist_ok=True)

#     plots = [
#         ("mc_plot1_idlmc_xrtpy_base.png",
#          lambda: plot1_idlmc_xrtpy_base(idl, xrtpy.dem, xrtpy.logT, mc_case)),
#         ("mc_plot2_idl_base_xrtpy_mc.png",
#          lambda: plot2_idl_base_xrtpy_mc(idl_true_base, idl.logT, xrtpy, mc_case)),
#         ("mc_plot3_both_ensembles.png",
#          lambda: plot3_both_ensembles(idl, xrtpy, mc_case, idl_true_base=idl_true_base)),
#         ("mc_plot4_combined.png",
#          lambda: plot4_combined(idl, xrtpy, mc_case, idl_true_base=idl_true_base)),
#     ]

#     for fname, make_fig in plots:
#         fig = make_fig()
#         if save:
#             out = case_plot_dir / fname
#             fig.savefig(out, dpi=180, bbox_inches="tight")
#             print(f"  Saved: {out}")
#             plt.close(fig)
#         else:
#             plt.show()
#             plt.close(fig)

#     print(f"{'='*60}")

# def main() -> None:
#     parser = argparse.ArgumentParser(
#         description="Generate MC DEM comparison plots (IDL vs XRTpy)."
#     )
#     parser.add_argument("sav_files", nargs="*",
#                         help="Specific MC .sav file(s) to plot.")
#     parser.add_argument("--all", action="store_true",
#                         help=f"Process all MC .sav files in {MC_DIR}")
#     parser.add_argument("--save", action="store_true",
#                         help="Save PNGs instead of showing interactively.")
#     args = parser.parse_args()

#     mc_cases: list[SavCase] = []

#     if args.all:
#         mc_cases = discover_mc_cases(MC_DIR)
#         if not mc_cases:
#             print(f"No MC .sav files found in {MC_DIR}")
#             sys.exit(1)
#         print(f"Found {len(mc_cases)} MC case(s)")
#     elif args.sav_files:
#         for p in args.sav_files:
#             mc_cases.append(parse_sav_filename(Path(p)))
#     else:
#         parser.print_help()
#         sys.exit(0)

#     for case in mc_cases:
#         try:
#             process_case(case, save=args.save)
#         except Exception as exc:
#             print(f"  ERROR processing {case.label}: {exc}")
#             raise

#     print("\nDone.")


# if __name__ == "__main__":
#     main()

"""
plot_mc_comparison.py
=====================
Four-panel Monte Carlo comparison: IDL vs XRTpy DEM ensembles.
Now includes a chi-square histogram panel (Plot 5 / Panel 4 in combined figure).

Usage
-----
    python plot_mc_comparison.py --case 20071213T0401
    python plot_mc_comparison.py --all
    python plot_mc_comparison.py --all --save

The script expects:
    data/validation/base/        xrt_IDL_dem_<LABEL>_..._.sav
    data/validation/monte_carlo/ xrt_IDL_dem_<LABEL>_..._MC<N>_.sav

Plots produced per case:
    plots/<LABEL>/mc_plot1_idlmc_xrtpy_base.png
    plots/<LABEL>/mc_plot2_idl_base_xrtpy_mc.png
    plots/<LABEL>/mc_plot3_both_ensembles.png
    plots/<LABEL>/mc_plot4_chisq_histogram.png
    plots/<LABEL>/mc_plot5_combined.png        <- shareable 4-panel figure
"""

import argparse
import sys
import warnings
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np
from utils_sav_io import (
    IDLMCResult,
    SavCase,
    discover_mc_cases,
    load_idl_mc_sav,
    load_idl_sav,
    parse_sav_filename,
)

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
BASE_DIR = Path(__file__).parent / "data" / "validation" / "base"
MC_DIR = Path(__file__).parent / "data" / "validation" / "monte_carlo"
PLOT_DIR = Path(__file__).parent / "plots"

IDL_COLOR = "#BA7517"  # amber
XRTPY_COLOR = "#185FA5"  # blue
IDL_ALPHA = 0.12
XRTPY_ALPHA = 0.12
BASE_LW = 2.5
ENS_LW = 0.6
DEM_FLOOR = 1e-99


def _log10(dem: np.ndarray) -> np.ndarray:
    return np.log10(np.maximum(dem, DEM_FLOOR))


def _axis_style(ax, title: str, xlim=(5.4, 8.1)):
    ax.set_xlabel("log₁₀ T  [K]", fontsize=10)
    ax.set_ylabel("log₁₀ DEM  [cm⁻⁵ K⁻¹]", fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.set_xlim(*xlim)
    ax.grid(True, alpha=0.22)


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------


def run_xrtpy(case: SavCase, n_mc: int) -> XRTDEMIterative:
    """Run XRTpy solver with the same inputs as the IDL case."""
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        responses = generate_temperature_responses(case.filters, case.observation_date)
        solver = XRTDEMIterative(
            observed_channel=case.filters,
            observed_intensities=case.intensities_array,
            temperature_responses=responses,
            monte_carlo_runs=n_mc,
        )
        solver.solve()
    return solver


# ---------------------------------------------------------------------------
# DEM Plot functions (1-3)
# ---------------------------------------------------------------------------


def plot1_idlmc_xrtpy_base(
    idl: IDLMCResult,
    xrtpy_base: np.ndarray,
    logT: np.ndarray,
    case: SavCase,
    idl_true_base: np.ndarray | None = None,
    ax=None,
) -> plt.Figure | None:
    """Plot 1: IDL 1000 MC ensemble + XRTpy base DEM."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 5.5))

    _base = idl_true_base if idl_true_base is not None else idl.dem_base

    for i in range(idl.n_runs):
        ax.step(
            idl.logT,
            _log10(idl.dem_mc[i]),
            where="mid",
            color=IDL_COLOR,
            alpha=IDL_ALPHA,
            linewidth=ENS_LW,
        )

    ax.step(
        idl.logT,
        _log10(_base),
        where="mid",
        color=IDL_COLOR,
        linewidth=BASE_LW,
        label="IDL base DEM",
    )
    ax.step(
        logT,
        _log10(xrtpy_base),
        where="mid",
        color=XRTPY_COLOR,
        linewidth=BASE_LW,
        linestyle="--",
        label="XRTpy base DEM",
    )

    _axis_style(
        ax,
        f"Plot 1 — IDL {idl.n_runs} MC + XRTpy base\n"
        f"Does XRTpy fall within IDL's MC spread?",
    )
    ax.legend(fontsize=9)

    if standalone:
        fig.suptitle(
            f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
            fontsize=10,
            y=1.01,
        )
        fig.tight_layout()
        return fig
    return None


def plot2_idl_base_xrtpy_mc(
    idl_base: np.ndarray,
    idl_logT: np.ndarray,
    xrtpy: XRTDEMIterative,
    case: SavCase,
    ax=None,
) -> plt.Figure | None:
    """Plot 2: IDL base DEM + XRTpy 1000 MC ensemble."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 5.5))

    n_mc = xrtpy.mc_dem.shape[0] - 1

    for i in range(1, n_mc + 1):
        ax.step(
            xrtpy.logT,
            _log10(xrtpy.mc_dem[i]),
            where="mid",
            color=XRTPY_COLOR,
            alpha=XRTPY_ALPHA,
            linewidth=ENS_LW,
        )

    ax.step(
        xrtpy.logT,
        _log10(xrtpy.dem),
        where="mid",
        color=XRTPY_COLOR,
        linewidth=BASE_LW,
        linestyle="--",
        label="XRTpy base DEM",
    )
    ax.step(
        idl_logT,
        _log10(idl_base),
        where="mid",
        color=IDL_COLOR,
        linewidth=BASE_LW,
        label="IDL base DEM",
    )

    _axis_style(
        ax,
        f"Plot 2 — IDL base + XRTpy {n_mc} MC\nDoes IDL fall within XRTpy's MC spread?",
    )
    ax.legend(fontsize=9)

    if standalone:
        fig.suptitle(
            f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
            fontsize=10,
            y=1.01,
        )
        fig.tight_layout()
        return fig
    return None


def plot3_both_ensembles(
    idl: IDLMCResult,
    xrtpy: XRTDEMIterative,
    case: SavCase,
    idl_true_base: np.ndarray | None = None,
    ax=None,
) -> plt.Figure | None:
    """Plot 3: IDL 1000 MC + XRTpy 1000 MC both overlaid."""
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 5.5))

    _base = idl_true_base if idl_true_base is not None else idl.dem_base
    n_xrt = xrtpy.mc_dem.shape[0] - 1

    for i in range(idl.n_runs):
        ax.step(
            idl.logT,
            _log10(idl.dem_mc[i]),
            where="mid",
            color=IDL_COLOR,
            alpha=IDL_ALPHA,
            linewidth=ENS_LW,
        )

    for i in range(1, n_xrt + 1):
        ax.step(
            xrtpy.logT,
            _log10(xrtpy.mc_dem[i]),
            where="mid",
            color=XRTPY_COLOR,
            alpha=XRTPY_ALPHA,
            linewidth=ENS_LW,
        )

    ax.step(
        idl.logT,
        _log10(_base),
        where="mid",
        color=IDL_COLOR,
        linewidth=BASE_LW,
        label=f"IDL base  ({idl.n_runs} MC shown)",
    )
    ax.step(
        xrtpy.logT,
        _log10(xrtpy.dem),
        where="mid",
        color=XRTPY_COLOR,
        linewidth=BASE_LW,
        linestyle="--",
        label=f"XRTpy base ({n_xrt} MC shown)",
    )

    _axis_style(
        ax,
        f"Plot 3 — Both ensembles overlaid\n"
        f"IDL {idl.n_runs} MC (amber) vs XRTpy {n_xrt} MC (blue)",
    )
    ax.legend(fontsize=9)

    if standalone:
        fig.suptitle(
            f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
            fontsize=10,
            y=1.01,
        )
        fig.tight_layout()
        return fig
    return None


# ---------------------------------------------------------------------------
# Chi-square histogram (Plot 4 / Panel 4)
# ---------------------------------------------------------------------------


def plot4_chisq_histogram(
    idl: IDLMCResult,
    xrtpy: XRTDEMIterative,
    case: SavCase,
    ax=None,
) -> plt.Figure | None:
    """
    Plot 4: Overlapping chi-square histograms — IDL (amber) vs XRTpy (blue).

    Shows:
    - Distribution of 1000 MC chi-square values for each code
    - Vertical dashed lines for base chi-square of each code
    - Mean annotations
    - Log x-axis (chi-square values span large ranges)

    Answers: Do IDL and XRTpy find similarly good fits across MC realizations?
    """
    standalone = ax is None
    if standalone:
        fig, ax = plt.subplots(figsize=(9, 5.0))

    # --- Gather chi-square arrays ---
    # XRTpy: row 0 = base, rows 1..N = MC
    xrt_chisq_mc = xrtpy.mc_chisq[1:].astype(float)
    xrt_chisq_base = float(xrtpy.mc_chisq[0])

    # IDL: use chisq_mc if available, otherwise skip IDL histogram
    idl_has_chisq = idl.chisq_mc is not None and len(idl.chisq_mc) > 0
    idl_chisq_mc = idl.chisq_mc.astype(float) if idl_has_chisq else None
    idl_chisq_base = float(idl.chisq_base)

    # --- Determine bin edges (shared log-spaced bins) ---
    all_vals = list(xrt_chisq_mc[xrt_chisq_mc > 0])
    if idl_has_chisq:
        all_vals += list(idl_chisq_mc[idl_chisq_mc > 0])

    if not all_vals:
        ax.text(
            0.5,
            0.5,
            "No chi-square data available",
            ha="center",
            va="center",
            transform=ax.transAxes,
        )
        if standalone:
            return plt.gcf()
        return None

    vmin = max(min(all_vals) * 0.5, 1.0)
    vmax = max(all_vals) * 2.0
    bins = np.logspace(np.log10(vmin), np.log10(vmax), 50)

    # --- Plot histograms ---
    if idl_has_chisq:
        ax.hist(
            idl_chisq_mc,
            bins=bins,
            color=IDL_COLOR,
            alpha=0.55,
            label=f"IDL MC (n={idl.n_runs})",
            edgecolor="none",
        )

    ax.hist(
        xrt_chisq_mc,
        bins=bins,
        color=XRTPY_COLOR,
        alpha=0.55,
        label=f"XRTpy MC (n={len(xrt_chisq_mc)})",
        edgecolor="none",
    )

    # --- Vertical lines for base chi-square ---
    if idl_has_chisq and idl_chisq_base > 0:
        ax.axvline(
            idl_chisq_base,
            color=IDL_COLOR,
            linewidth=2.0,
            linestyle="--",
            label=f"IDL base χ² = {idl_chisq_base:.1f}",
        )

    ax.axvline(
        xrt_chisq_base,
        color=XRTPY_COLOR,
        linewidth=2.0,
        linestyle="--",
        label=f"XRTpy base χ² = {xrt_chisq_base:.1f}",
    )

    # --- Annotate means ---
    xrt_mean = float(np.mean(xrt_chisq_mc))
    ax.axvline(xrt_mean, color=XRTPY_COLOR, linewidth=1.2, linestyle=":", alpha=0.8)
    ax.text(
        xrt_mean * 1.05,
        ax.get_ylim()[1] * 0.85 if ax.get_ylim()[1] > 0 else 1,
        f"mean={xrt_mean:.1f}",
        color=XRTPY_COLOR,
        fontsize=8,
        va="top",
    )

    if idl_has_chisq:
        idl_mean = float(np.mean(idl_chisq_mc))
        ax.axvline(idl_mean, color=IDL_COLOR, linewidth=1.2, linestyle=":", alpha=0.8)
        ax.text(
            idl_mean * 1.05,
            ax.get_ylim()[1] * 0.70 if ax.get_ylim()[1] > 0 else 1,
            f"mean={idl_mean:.1f}",
            color=IDL_COLOR,
            fontsize=8,
            va="top",
        )

    # --- Axis formatting ---
    ax.set_xscale("log")
    ax.set_xlabel("χ²", fontsize=10)
    ax.set_ylabel("Count", fontsize=10)
    ax.set_title(
        "Plot 4 — χ² distribution: IDL vs XRTpy\n"
        "Are the fits equally good across MC realizations?",
        fontsize=11,
    )
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(True, alpha=0.22, which="both")

    # Annotate with summary stats in a text box
    if idl_has_chisq:
        summary = (
            f"IDL:   median={np.median(idl_chisq_mc):.1f}  σ={np.std(idl_chisq_mc):.1f}\n"
            f"XRTpy: median={np.median(xrt_chisq_mc):.1f}  σ={np.std(xrt_chisq_mc):.1f}"
        )
    else:
        summary = (
            f"XRTpy: median={np.median(xrt_chisq_mc):.1f}  σ={np.std(xrt_chisq_mc):.1f}\n"
            f"IDL chi-square not available in SAV"
        )
    ax.text(
        0.02,
        0.97,
        summary,
        transform=ax.transAxes,
        va="top",
        ha="left",
        fontsize=8.5,
        bbox=dict(
            boxstyle="round,pad=0.4", facecolor="white", edgecolor="gray", alpha=0.85
        ),
    )

    if standalone:
        fig.suptitle(
            f"Case: {case.observation_date}  |  Filters: {', '.join(case.filters)}",
            fontsize=10,
            y=1.01,
        )
        fig.tight_layout()
        return fig
    return None


# ---------------------------------------------------------------------------
# Combined 4-panel figure (Plot 5)
# ---------------------------------------------------------------------------


def plot5_combined(
    idl: IDLMCResult,
    xrtpy: XRTDEMIterative,
    case: SavCase,
    idl_true_base: np.ndarray | None = None,
) -> plt.Figure:
    """
    Plot 5: All four panels stacked — the shareable team figure.

    Panel layout:
        1. IDL MC + XRTpy base
        2. IDL base + XRTpy MC
        3. Both ensembles overlaid
        4. Chi-square histogram
    """
    fig = plt.figure(figsize=(10, 22))

    fig.suptitle(
        f"IDL vs XRTpy — Monte Carlo DEM Comparison\n"
        f"Case: {case.observation_date}   "
        f"Filters: {', '.join(case.filters)}\n"
        f"IDL: {idl.n_runs} MC runs   |   "
        f"XRTpy: {xrtpy.mc_dem.shape[0] - 1} MC runs",
        fontsize=12,
        fontweight="bold",
        y=0.99,
    )

    gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.46)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2])
    ax4 = fig.add_subplot(gs[3])

    _base = idl_true_base if idl_true_base is not None else idl.dem_base

    plot1_idlmc_xrtpy_base(
        idl, xrtpy.dem, xrtpy.logT, case, idl_true_base=_base, ax=ax1
    )
    plot2_idl_base_xrtpy_mc(_base, idl.logT, xrtpy, case, ax=ax2)
    plot3_both_ensembles(idl, xrtpy, case, idl_true_base=_base, ax=ax3)
    plot4_chisq_histogram(idl, xrtpy, case, ax=ax4)

    return fig


# ---------------------------------------------------------------------------
# Main process function
# ---------------------------------------------------------------------------


def process_case(mc_case: SavCase, save: bool) -> None:
    """Load data, run XRTpy, generate all 5 plots for one case."""
    print(f"\n{'=' * 60}")
    print(f"  Case: {mc_case.label}  ({mc_case.observation_date})")
    print(f"  Filters: {mc_case.filters}")
    print(f"  IDL MC runs: {mc_case.mc_runs}")

    # Load IDL MC results
    print(f"  Loading MC SAV: {mc_case.sav_path.name}")
    idl = load_idl_mc_sav(mc_case.sav_path)
    print(f"  IDL dem_mc shape: {idl.dem_mc.shape}  (n_runs x nT)")

    # Chi-square status
    if idl.chisq_mc is not None:
        print(
            f"  IDL chisq_mc: shape={idl.chisq_mc.shape}  "
            f"median={np.median(idl.chisq_mc):.1f}  base={idl.chisq_base:.1f}"
        )
    else:
        print(f"  IDL chisq_mc: not available in SAV  base={idl.chisq_base:.1f}")

    # Load true base DEM from standalone base SAV
    base_matches = [
        f
        for f in BASE_DIR.glob(f"xrt_IDL_dem_{mc_case.label}_*.sav")
        if "MC" not in f.name
    ]
    if base_matches:
        idl_base_result = load_idl_sav(base_matches[0])
        idl_true_base = idl_base_result.dem
        print(f"  Loaded true base DEM from: {base_matches[0].name}")
    else:
        print(f"  WARNING: no base SAV found for {mc_case.label} — using MC row 0")
        idl_true_base = idl.dem_base

    # Run XRTpy
    print(f"  Running XRTpy with {mc_case.mc_runs} MC runs...")
    xrtpy = run_xrtpy(mc_case, n_mc=mc_case.mc_runs)
    print(f"  XRTpy mc_dem shape: {xrtpy.mc_dem.shape}")
    print(f"  XRTpy base χ²: {xrtpy.chisq:.2f}")
    print(f"  XRTpy MC χ² median: {np.median(xrtpy.mc_chisq[1:]):.2f}")

    # Output folder
    case_plot_dir = PLOT_DIR / mc_case.label
    case_plot_dir.mkdir(parents=True, exist_ok=True)

    plots = [
        (
            "mc_plot1_idlmc_xrtpy_base.png",
            lambda: plot1_idlmc_xrtpy_base(
                idl, xrtpy.dem, xrtpy.logT, mc_case, idl_true_base=idl_true_base
            ),
        ),
        (
            "mc_plot2_idl_base_xrtpy_mc.png",
            lambda: plot2_idl_base_xrtpy_mc(idl_true_base, idl.logT, xrtpy, mc_case),
        ),
        (
            "mc_plot3_both_ensembles.png",
            lambda: plot3_both_ensembles(
                idl, xrtpy, mc_case, idl_true_base=idl_true_base
            ),
        ),
        (
            "mc_plot4_chisq_histogram.png",
            lambda: plot4_chisq_histogram(idl, xrtpy, mc_case),
        ),
        (
            "mc_plot5_combined.png",
            lambda: plot5_combined(idl, xrtpy, mc_case, idl_true_base=idl_true_base),
        ),
    ]

    for fname, make_fig in plots:
        fig = make_fig()
        if save:
            out = case_plot_dir / fname
            fig.savefig(out, dpi=180, bbox_inches="tight")
            print(f"  Saved: {out}")
            plt.close(fig)
        else:
            plt.show()
            plt.close(fig)

    print(f"{'=' * 60}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate MC DEM comparison plots (IDL vs XRTpy)."
    )
    parser.add_argument(
        "sav_files", nargs="*", help="Specific MC .sav file(s) to plot."
    )
    parser.add_argument(
        "--all", action="store_true", help=f"Process all MC .sav files in {MC_DIR}"
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save PNGs instead of showing interactively.",
    )
    args = parser.parse_args()

    mc_cases: list[SavCase] = []

    if args.all:
        mc_cases = discover_mc_cases(MC_DIR)
        if not mc_cases:
            print(f"No MC .sav files found in {MC_DIR}")
            sys.exit(1)
        print(f"Found {len(mc_cases)} MC case(s)")
    elif args.sav_files:
        for p in args.sav_files:
            mc_cases.append(parse_sav_filename(Path(p)))
    else:
        parser.print_help()
        sys.exit(0)

    for case in mc_cases:
        try:
            process_case(case, save=args.save)
        except Exception as exc:
            print(f"  ERROR processing {case.label}: {exc}")
            raise

    print("\nDone.")


if __name__ == "__main__":
    main()

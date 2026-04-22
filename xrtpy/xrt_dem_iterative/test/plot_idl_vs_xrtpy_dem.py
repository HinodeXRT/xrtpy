"""
IDL vs XRTpy DEM Comparison — Multi-Case Plot Script
======================================================

Usage
-----
Plot one specific .sav file:
    python plot_idl_vs_xrtpy_dem.py data/validation/xrt_IDL_dem_20080104T1104_..._.sav

Plot all .sav files in the validation directory:
    python plot_idl_vs_xrtpy_dem.py --all

Plot all and save as PNG (no interactive window):
    python plot_idl_vs_xrtpy_dem.py --all --save

The script reads filter names and intensities directly from the filename —
no hardcoding required.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from xrtpy.xrt_dem_iterative.utils_sav_io import SavCase, IDLResult, discover_cases, load_idl_sav, parse_sav_filename

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
MEAN_TOL  = 0.20
MAX_TOL   = 0.50
DEM_FLOOR = 1e10


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _log10_safe(arr: np.ndarray, floor: float = 1e-99) -> np.ndarray:
    return np.log10(np.maximum(arr, floor))


def _valid_mask(dem_idl: np.ndarray, dem_xrt: np.ndarray) -> np.ndarray:
    return (dem_idl > DEM_FLOOR) | (dem_xrt > DEM_FLOOR)


def _bar_color(d: float) -> str:
    a = abs(d)
    if a <= MEAN_TOL:
        return "#1D9E75"
    if a <= MAX_TOL:
        return "#BA7517"
    return "#E24B4A"


# ---------------------------------------------------------------------------
# Core: solve one case
# ---------------------------------------------------------------------------

def solve_case(case: SavCase) -> tuple[IDLResult, XRTDEMIterative]:
    print(f"  Loading IDL .sav: {case.sav_path.name}")
    idl = load_idl_sav(case.sav_path)

    print(f"  Running XRTpy solver ({len(case.filters)} filters)...")
    responses = generate_temperature_responses(case.filters, case.observation_date)
    solver = XRTDEMIterative(
        observed_channel=case.filters,
        observed_intensities=case.intensities_array,
        temperature_responses=responses,
        monte_carlo_runs=0,
    )
    solver.solve()

    return idl, solver


# ---------------------------------------------------------------------------
# Core: make the 3-panel comparison figure
# ---------------------------------------------------------------------------

def make_figure(case: SavCase, idl: IDLResult, solver: XRTDEMIterative) -> plt.Figure:

    logT_idl = idl.logT
    log_dem_idl = _log10_safe(idl.dem)

    logT_xrt = solver.logT
    log_dem_xrt = _log10_safe(solver.dem)
    modeled = solver.modeled_intensities
    chisq = solver.chisq

    mask = _valid_mask(idl.dem, solver.dem)
    delta = log_dem_xrt - log_dem_idl
    mean_diff = float(np.mean(np.abs(delta[mask])))
    max_diff  = float(np.max(np.abs(delta[mask])))

    log_ratio = np.log10(np.maximum(modeled / case.intensities_array, 1e-99))

    # Print summary to terminal
    print(f"\n  {'─'*55}")
    print(f"  Case:          {case.label}  ({case.observation_date})")
    print(f"  Filters:       {case.filters}")
    print(f"  XRTpy χ²:      {chisq:.2f}")
    print(f"  Mean |Δ|:      {mean_diff:.4f} dex")
    print(f"  Max  |Δ|:      {max_diff:.4f} dex")
    print(f"  IDL  peak:     logT = {logT_idl[np.argmax(idl.dem)]:.2f}")
    print(f"  XRTpy peak:    logT = {logT_xrt[np.argmax(solver.dem)]:.2f}")
    print(f"  Per-filter log10(modeled/observed):")
    for f, lr in zip(case.filters, log_ratio):
        flag = "  ← !" if abs(lr) > 1.0 else ""
        print(f"    {f:<22} {lr:+.3f}{flag}")
    print(f"  {'─'*55}")

    # ── Build figure ─────────────────────────────────────────────────────
    fig = plt.figure(figsize=(11, 13))
    fig.suptitle(
        f"IDL vs XRTpy DEM Comparison\n"
        f"Case: {case.observation_date}   Filters: {', '.join(case.filters)}",
        fontsize=12,
        fontweight="bold",
        y=0.98,
    )

    gs = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[3, 2, 2], hspace=0.48)

    # ── Panel 1: DEM curves ──────────────────────────────────────────────
    ax1 = fig.add_subplot(gs[0])

    ax1.step(logT_idl, log_dem_idl, where="mid",
             color="#BA7517", linewidth=2.2, label="IDL reference")
    ax1.step(logT_xrt, log_dem_xrt, where="mid",
             color="#185FA5", linewidth=2.0, linestyle="--", label="XRTpy")

    # Shade high-T region where secondary rise often lives
    ax1.axvspan(7.3, 8.05, alpha=0.07, color="crimson", label="High-T region (>logT 7.3)")

    # Mark spline knot positions
    knot_logT = np.linspace(logT_xrt.min(), logT_xrt.max(), solver.n_spl)
    for j, klt in enumerate(knot_logT):
        ax1.axvline(
            klt,
            color="#3B8BD4",
            linewidth=0.7,
            alpha=0.5,
            linestyle=":",
            label=f"XRTpy spline knots (n={solver.n_spl})" if j == 0 else "_nolegend_",
        )

    ax1.set_xlabel("log₁₀ T  [K]")
    ax1.set_ylabel("log₁₀ DEM  [cm⁻⁵ K⁻¹]")
    ax1.set_title("DEM(T) comparison", fontsize=11)
    ax1.set_xlim(5.4, 8.1)
    ax1.legend(fontsize=9, loc="upper right")
    ax1.grid(True, alpha=0.25)

    # Annotate peaks
    pk_idl = logT_idl[np.argmax(idl.dem)]
    pk_xrt = logT_xrt[np.argmax(solver.dem)]
    y_top_idl = log_dem_idl.max()
    y_top_xrt = log_dem_xrt.max()

    #Maybelater 
    # ax1.annotate(
    #     f"IDL\nlogT={pk_idl:.2f}",
    #     xy=(pk_idl, y_top_idl),
    #     xytext=(pk_idl - 0.7, y_top_idl - 2.5),
    #     fontsize=8, color="#BA7517",
    #     arrowprops=dict(arrowstyle="->", color="#BA7517", lw=0.8),
    # )
    # ax1.annotate(
    #     f"XRTpy\nlogT={pk_xrt:.2f}",
    #     xy=(pk_xrt, y_top_xrt),
    #     xytext=(pk_xrt + 0.15, y_top_xrt - 3.5),
    #     fontsize=8, color="#185FA5",
    #     arrowprops=dict(arrowstyle="->", color="#185FA5", lw=0.8),
    # )

    # ── Panel 2: Δ per bin ───────────────────────────────────────────────
    ax2 = fig.add_subplot(gs[1])

    bar_colors = [_bar_color(d) for d in delta]
    ax2.bar(logT_xrt, delta, width=0.09, color=bar_colors, alpha=0.85)

    ax2.axhline(0, color="black", linewidth=1.0, alpha=0.35)
    ax2.axhline(+MEAN_TOL, color="#1D9E75", linewidth=1.0, linestyle="--",
                label=f"±{MEAN_TOL} dex (mean tol.)", alpha=0.8)
    ax2.axhline(-MEAN_TOL, color="#1D9E75", linewidth=1.0, linestyle="--", alpha=0.8)
    ax2.axhline(+MAX_TOL,  color="#BA7517", linewidth=1.0, linestyle=":",
                label=f"±{MAX_TOL} dex (max tol.)", alpha=0.8)
    ax2.axhline(-MAX_TOL,  color="#BA7517", linewidth=1.0, linestyle=":", alpha=0.8)
    ax2.fill_between([5.4, 8.1], -MEAN_TOL, MEAN_TOL, alpha=0.06, color="#1D9E75")

    ax2.set_xlabel("log₁₀ T  [K]")
    ax2.set_ylabel("Δ log₁₀ DEM\n(XRTpy − IDL)")
    ax2.set_xlim(5.4, 8.1)
    ax2.set_title(
        f"Per-bin difference   mean|Δ|={mean_diff:.3f} dex   max|Δ|={max_diff:.3f} dex",
        fontsize=11,
    )
    ax2.legend(fontsize=9, loc="upper left")
    ax2.grid(True, alpha=0.25, axis="y")

    # ── Panel 3: per-filter intensities ──────────────────────────────────
    ax3 = fig.add_subplot(gs[2])

    n_filters = len(case.filters)
    x = np.arange(n_filters)
    w = 0.35

    ax3.bar(x - w / 2, case.intensities_array, w,
            label="Observed", color="#185FA5", alpha=0.85)
    ax3.bar(x + w / 2, modeled, w,
            label="XRTpy modeled", color="#85B7EB", alpha=0.85)

    for i, (obs, mod) in enumerate(zip(case.intensities_array, modeled)):
        lr = np.log10(max(mod / obs, 1e-99))
        color = "#E24B4A" if abs(lr) > 1.0 else "#185FA5"
        ax3.text(
            i + w / 2,
            mod * 1.3,
            f"{lr:+.2f}",
            ha="center",
            va="bottom",
            fontsize=8.5,
            color=color,
            fontweight="bold",
        )

    ax3.set_yscale("log")
    ax3.set_xticks(x)
    ax3.set_xticklabels(case.filters, rotation=20, ha="right", fontsize=9)
    ax3.set_ylabel("Intensity  [DN/s/pix]")
    ax3.set_title(
        "Per-filter intensities — numbers show log₁₀(modeled / observed)",
        fontsize=11,
    )
    ax3.legend(fontsize=9)
    ax3.grid(True, alpha=0.25, axis="y")
    ax3.set_xlim(-0.6, n_filters - 0.4)

    ax3.text(
        0.98, 0.97,
        f"XRTpy χ² = {chisq:.1f}",
        transform=ax3.transAxes,
        ha="right", va="top",
        fontsize=9, color="#E24B4A",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                  edgecolor="#E24B4A", alpha=0.8),
    )

    return fig


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description="Plot IDL vs XRTpy DEM comparison for one or more .sav files."
    )
    parser.add_argument(
        "sav_files",
        nargs="*",
        help="Path(s) to specific .sav file(s) to plot.",
    )
    parser.add_argument(
        "--all",
        action="store_true",
        help="Plot all xrt_IDL_dem_*.sav files in data/validation/.",
    )
    parser.add_argument(
        "--save",
        action="store_true",
        help="Save each plot as a PNG file instead of showing interactively.",
    )
    args = parser.parse_args()

    # Collect cases to plot
    cases: list[SavCase] = []

    if args.all:
        data_dir = Path(__file__).parent / "data" / "validation"
        cases = discover_cases(data_dir)
        if not cases:
            print(f"No xrt_IDL_dem_*.sav files found in {data_dir}")
            sys.exit(1)
        print(f"Found {len(cases)} case(s) in {data_dir}")
    elif args.sav_files:
        for p in args.sav_files:
            cases.append(parse_sav_filename(Path(p)))
    else:
        parser.print_help()
        sys.exit(0)

    # Process each case
    for i, case in enumerate(cases, 1):
        print(f"\n[{i}/{len(cases)}] Case: {case.label}")
        print(f"  Date:    {case.observation_date}")
        print(f"  Filters: {case.filters}")
        print(f"  I_obs:   {[f'{v:.3f}' for v in case.intensities]}")

        try:
            idl, solver = solve_case(case)
        except Exception as exc:
            print(f"  ERROR: {exc}")
            continue

        fig = make_figure(case, idl, solver)

        # if args.save:
        #     out_path = Path(f"dem_comparison_{case.label}.png")
        #     fig.savefig(out_path, dpi=180, bbox_inches="tight")
        #     print(f"  Saved: {out_path}")
        #     plt.close(fig)
        if args.save:
            case_dir = Path("plots") / case.label
            case_dir.mkdir(parents=True, exist_ok=True)
            out_path = case_dir / "base_dem.png"
            fig.savefig(out_path, dpi=180, bbox_inches="tight")
            print(f"  Saved: {out_path}")
            plt.close(fig)
        else:
            plt.tight_layout()
            plt.show()

    print("\nDone.")


if __name__ == "__main__":
    main()

# """
# IDL vs XRTpy DEM Comparison Plot
# =================================
# Run this script to generate the real comparison using your actual XRTpy solver
# and the IDL reference .sav file.

# Usage:
#     python plot_idl_vs_xrtpy_dem.py

# Requires:
#     - scipy, numpy, matplotlib
#     - xrtpy installed
#     - xrt_IDL_dem_output_20080104.sav in the same directory
# """

# from pathlib import Path

# import matplotlib.gridspec as gridspec
# import matplotlib.pyplot as plt
# import numpy as np
# from scipy.io import readsav

# from xrtpy.response.tools import generate_temperature_responses
# from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # ---------------------------------------------------------------------------
# # Config
# # ---------------------------------------------------------------------------
# DATA_DIR = Path(__file__).parent / "data" / "validation"

# SAV_FILE = DATA_DIR / "xrt_dem_output_20071213T0401_NOMC.sav"#"xrt_IDL_dem_output_20080104.sav"

# # FILTERS = ["Be-med", "Al-mesh", "Ti-poly", "Al-poly", "Be-thin"]
# # INTENSITIES = np.array([234.283365, 183.711876, 45.931438, 91.745329, 5.755926])
# # OBSERVATION_DATE = "2008-01-04T11:04:26"
    
# FILTERS= ["be-med","Be-thin","Al-poly", "Al-poly/Ti-poly","Ti-poly","Al-thick"]
# INTENSITIES = [603.875886,150.921435,2412.340960, 301.354389 ,603.100596,2.519851]
# OBSERVATION_DATE = "2007-12-13T04:01"

    
# DEM_FLOOR = 1e10

# MEAN_TOL = 0.20
# MAX_TOL  = 0.50

# # ---------------------------------------------------------------------------
# # Load IDL reference
# # ---------------------------------------------------------------------------
# print("Loading IDL .sav file...")
# data = readsav(str(SAV_FILE), python_dict=True)
# logT_idl = np.array(data["logt"]).ravel().astype(float)
# dem_idl  = np.array(data["dem"]).ravel().astype(float)
# log_dem_idl = np.log10(np.maximum(dem_idl, 1e-99))

# # ---------------------------------------------------------------------------
# # Run XRTpy solver
# # ---------------------------------------------------------------------------
# print("Running XRTpy DEM solver...")
# responses = generate_temperature_responses(FILTERS, OBSERVATION_DATE)
# solver = XRTDEMIterative(
#     observed_channel=FILTERS,
#     observed_intensities=INTENSITIES,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# solver.solve()

# logT_xrt    = solver.logT
# dem_xrt     = solver.dem
# log_dem_xrt = np.log10(np.maximum(dem_xrt, 1e-99))
# modeled     = solver.modeled_intensities
# chisq       = solver.chisq

# print(f"  XRTpy chi-square = {chisq:.2f}")
# print(f"  IDL  peak logT   = {logT_idl[np.argmax(dem_idl)]:.2f}")
# print(f"  XRTpy peak logT  = {logT_xrt[np.argmax(dem_xrt)]:.2f}")

# # ---------------------------------------------------------------------------
# # Metrics
# # ---------------------------------------------------------------------------
# mask = (dem_idl > DEM_FLOOR) | (dem_xrt > DEM_FLOOR)
# delta = log_dem_xrt - log_dem_idl
# mean_diff = np.mean(np.abs(delta[mask]))
# max_diff  = np.max(np.abs(delta[mask]))
# print(f"  Mean |Δlog10(DEM)| = {mean_diff:.3f} dex")
# print(f"  Max  |Δlog10(DEM)| = {max_diff:.3f} dex")

# log_ratio = np.log10(np.maximum(modeled / INTENSITIES, 1e-99))
# print(f"  Per-filter log10(modeled/observed):")
# for f, lr in zip(FILTERS, log_ratio):
#     flag = "  ← !" if abs(lr) > 1.0 else ""
#     print(f"    {f:<18} {lr:+.3f} dex{flag}")

# # ---------------------------------------------------------------------------
# # Plot
# # ---------------------------------------------------------------------------
# fig = plt.figure(figsize=(11, 13))
# fig.suptitle(
#     f"IDL vs XRTpy DEM Comparison\n"
#     f"Case: {OBSERVATION_DATE}   Filters: {', '.join(FILTERS)}",
#     fontsize=13, fontweight="bold", y=0.98
# )

# gs = gridspec.GridSpec(
#     3, 1,
#     figure=fig,
#     height_ratios=[3, 2, 2],
#     hspace=0.45
# )

# # ── Panel 1: DEM curves ──────────────────────────────────────────────────
# ax1 = fig.add_subplot(gs[0])

# ax1.step(logT_idl, log_dem_idl, where="mid",
#          color="#BA7517", linewidth=2.2, label="IDL reference")
# ax1.step(logT_xrt, log_dem_xrt, where="mid",
#          color="#185FA5", linewidth=2.0, linestyle="--", label="XRTpy")

# # Highlight the secondary rise region
# ax1.axvspan(7.3, 8.0, alpha=0.08, color="crimson", label="Missing high-T region")

# # Mark knot positions
# knot_logT = np.linspace(logT_xrt.min(), logT_xrt.max(), solver.n_spl)
# for klt in knot_logT:
#     ax1.axvline(klt, color="#3B8BD4", linewidth=0.7, alpha=0.5, linestyle=":")
# ax1.axvline(knot_logT[0], color="#3B8BD4", linewidth=0.7, alpha=0.5,
#             linestyle=":", label=f"XRTpy spline knots (n={solver.n_spl})")

# ax1.set_xlabel("log₁₀ T  [K]")
# ax1.set_ylabel("log₁₀ DEM  [cm⁻⁵ K⁻¹]")
# ax1.set_title("DEM(T) comparison", fontsize=11)
# ax1.set_xlim(5.4, 8.1)
# ax1.legend(fontsize=9, loc="upper right")
# ax1.grid(True, alpha=0.25)

# # Annotate peak
# pk_idl = logT_idl[np.argmax(dem_idl)]
# pk_xrt = logT_xrt[np.argmax(dem_xrt)]
# ax1.annotate(f"IDL peak\nlogT={pk_idl:.2f}",
#              xy=(pk_idl, log_dem_idl.max()),
#              xytext=(pk_idl - 0.6, log_dem_idl.max() - 3),
#              fontsize=8, color="#BA7517",
#              arrowprops=dict(arrowstyle="->", color="#BA7517", lw=0.8))
# ax1.annotate(f"XRTpy peak\nlogT={pk_xrt:.2f}",
#              xy=(pk_xrt, log_dem_xrt.max()),
#              xytext=(pk_xrt + 0.15, log_dem_xrt.max() - 3.5),
#              fontsize=8, color="#185FA5",
#              arrowprops=dict(arrowstyle="->", color="#185FA5", lw=0.8))

# # ── Panel 2: Delta per bin ───────────────────────────────────────────────
# ax2 = fig.add_subplot(gs[1])

# bar_colors = []
# for d in delta:
#     a = abs(d)
#     if a <= MEAN_TOL:
#         bar_colors.append("#1D9E75")
#     elif a <= MAX_TOL:
#         bar_colors.append("#BA7517")
#     else:
#         bar_colors.append("#E24B4A")

# ax2.bar(logT_xrt, delta, width=0.09, color=bar_colors, alpha=0.8)

# ax2.axhline(0,           color="black",   linewidth=1.0, alpha=0.4)
# ax2.axhline(+MEAN_TOL,   color="#1D9E75", linewidth=1.0, linestyle="--",
#             label=f"±{MEAN_TOL} dex (mean tol.)", alpha=0.8)
# ax2.axhline(-MEAN_TOL,   color="#1D9E75", linewidth=1.0, linestyle="--", alpha=0.8)
# ax2.axhline(+MAX_TOL,    color="#BA7517", linewidth=1.0, linestyle=":",
#             label=f"±{MAX_TOL} dex (max tol.)", alpha=0.8)
# ax2.axhline(-MAX_TOL,    color="#BA7517", linewidth=1.0, linestyle=":", alpha=0.8)

# ax2.fill_between([5.4, 8.1], -MEAN_TOL, MEAN_TOL, alpha=0.06, color="#1D9E75")

# ax2.set_xlabel("log₁₀ T  [K]")
# ax2.set_ylabel("Δ log₁₀ DEM\n(XRTpy − IDL)")
# ax2.set_xlim(5.4, 8.1)
# ax2.set_title(
#     f"Per-bin difference   mean|Δ|={mean_diff:.3f} dex   max|Δ|={max_diff:.3f} dex",
#     fontsize=11
# )
# ax2.legend(fontsize=9, loc="upper left")
# ax2.grid(True, alpha=0.25, axis="y")

# # ── Panel 3: Per-filter intensities ─────────────────────────────────────
# ax3 = fig.add_subplot(gs[2])

# x = np.arange(len(FILTERS))
# w = 0.35

# b1 = ax3.bar(x - w/2, INTENSITIES, w, label="Observed",
#              color="#185FA5", alpha=0.85)
# b2 = ax3.bar(x + w/2, modeled, w, label="XRTpy modeled",
#              color="#85B7EB", alpha=0.85)

# # Annotate log-ratio on each modeled bar
# for i, (obs, mod) in enumerate(zip(INTENSITIES, modeled)):
#     ratio = np.log10(max(mod / obs, 1e-99))
#     color = "#E24B4A" if abs(ratio) > 1.0 else "#3B8BD4"
#     ax3.text(
#         i + w/2, mod * 1.15,
#         f"{ratio:+.2f}",
#         ha="center", va="bottom",
#         fontsize=8.5, color=color, fontweight="bold"
#     )

# ax3.set_yscale("log")
# ax3.set_xticks(x)
# ax3.set_xticklabels(FILTERS)
# ax3.set_ylabel("Intensity  [DN/s/pix]")
# ax3.set_title(
#     "Per-filter intensities — numbers show log₁₀(modeled/observed)",
#     fontsize=11
# )
# ax3.legend(fontsize=9)
# ax3.grid(True, alpha=0.25, axis="y")
# ax3.set_xlim(-0.5, len(FILTERS) - 0.5)

# # Add chi-square annotation
# ax3.text(
#     0.98, 0.97,
#     f"XRTpy χ² = {chisq:.1f}",
#     transform=ax3.transAxes,
#     ha="right", va="top",
#     fontsize=9, color="#E24B4A",
#     bbox=dict(boxstyle="round,pad=0.3", facecolor="white", edgecolor="#E24B4A", alpha=0.8)
# )

# # ---------------------------------------------------------------------------
# # Save
# # ---------------------------------------------------------------------------
# out_path = Path("dem_idl_vs_xrtpy_comparison.png")
# fig.savefig(out_path, dpi=180, bbox_inches="tight")
# print(f"\nPlot saved to: {out_path}")
# plt.show()
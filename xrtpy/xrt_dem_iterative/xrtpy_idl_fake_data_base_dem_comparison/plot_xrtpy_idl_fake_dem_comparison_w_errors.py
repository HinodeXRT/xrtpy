# plot_base_dem_comparison_w_errors.py
#
# Plots three curves per run:
#   - Black solid   : True DEM (from DEM_XX.txt)
#   - Orange solid  : IDL base DEM (from IDL .sav with errors)
#   - Blue dashed   : XRTpy base DEM (solved with intensity_uncertainties)
#
# 24 runs total: 6 DEMs x 4 filter combos (4, 5, 6, 7 filters)
# Output: xrtpy_idl_fake_data_base_dem_comparison/

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR = Path("/Users/jvelasq/Projects/xrtpy/xrtpy/xrt_dem_iterative")
TRUE_DEM_DIR = BASE_DIR / "Fake_IDL_data" / "Fake_Based_IDL_DEMs"
IDL_SAV_DIR = (
    BASE_DIR / "xrtpy_idl_fake_data_base_dem_comparison" / "IDL_fake_DEM_data_w_errors"
)
OUT_DIR = (
    BASE_DIR
    / "xrtpy_idl_fake_data_base_dem_comparison"
    / "XRTpy_IDL_Fake_comparison_w_errors_overplots"
)
OUT_DIR.mkdir(exist_ok=True)
# OUT_DIR      = BASE_DIR / "xrtpy_idl_fake_data_base_dem_comparison"

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL = "10-Oct-2008 00:00:00"

# ── All 24 runs ───────────────────────────────────────────────────────────────
RUNS = [
    # ── DEM 14 ────────────────────────────────────────────────────────────────
    dict(
        dem_id=14,
        nfilters=4,
        filters=["Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[1112.14, 216.758, 94.3998, 8.71282],
        errors=[33.9430, 6.29030, 2.94368, 0.380063],
    ),
    dict(
        dem_id=14,
        nfilters=5,
        filters=["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[2329.27, 1112.14, 216.758, 94.3998, 8.71282],
        errors=[69.4233, 33.9430, 6.29030, 2.94368, 0.380063],
    ),
    dict(
        dem_id=14,
        nfilters=6,
        filters=["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[3584.44, 2329.27, 1112.14, 216.758, 94.3998, 8.71282],
        errors=[102.299, 69.4233, 33.9430, 6.29030, 2.94368, 0.380063],
    ),
    dict(
        dem_id=14,
        nfilters=7,
        filters=[
            "C-poly",
            "Ti-poly",
            "Be-thin",
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
        ],
        intensities=[3584.44, 2329.27, 1112.14, 216.758, 94.3998, 8.71282, 0.329589],
        errors=[102.299, 69.4233, 33.9430, 6.29030, 2.94368, 0.380063, 0.0873814],
    ),
    # ── DEM 15 ────────────────────────────────────────────────────────────────
    dict(
        dem_id=15,
        nfilters=4,
        filters=["Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[2633.15, 593.411, 243.889, 23.1174],
        errors=[73.6869, 17.5156, 6.66166, 0.728584],
    ),
    dict(
        dem_id=15,
        nfilters=5,
        filters=["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[4310.72, 2633.15, 593.411, 243.889, 23.1174],
        errors=[133.715, 73.6869, 17.5156, 6.66166, 0.728584],
    ),
    dict(
        dem_id=15,
        nfilters=6,
        filters=["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[6206.10, 4310.72, 2633.15, 593.411, 243.889, 23.1174],
        errors=[190.757, 133.715, 73.6869, 17.5156, 6.66166, 0.728584],
    ),
    dict(
        dem_id=15,
        nfilters=7,
        filters=[
            "C-poly",
            "Ti-poly",
            "Be-thin",
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
        ],
        intensities=[6206.10, 4310.72, 2633.15, 593.411, 243.889, 23.1174, 1.48818],
        errors=[190.757, 133.715, 73.6869, 17.5156, 6.66166, 0.728584, 0.166503],
    ),
    # ── DEM 65 ────────────────────────────────────────────────────────────────
    dict(
        dem_id=65,
        nfilters=4,
        filters=["Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[4331.36, 1007.78, 416.894, 40.2357],
        errors=[134.025, 27.1269, 12.3450, 1.13894],
    ),
    dict(
        dem_id=65,
        nfilters=5,
        filters=["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[7020.42, 4331.36, 1007.78, 416.894, 40.2357],
        errors=[202.527, 134.025, 27.1269, 12.3450, 1.13894],
    ),
    dict(
        dem_id=65,
        nfilters=6,
        filters=["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[10060.6, 7020.42, 4331.36, 1007.78, 416.894, 40.2357],
        errors=[288.248, 202.527, 134.025, 27.1269, 12.3450, 1.13894],
    ),
    dict(
        dem_id=65,
        nfilters=7,
        filters=[
            "C-poly",
            "Ti-poly",
            "Be-thin",
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
        ],
        intensities=[10060.6, 7020.42, 4331.36, 1007.78, 416.894, 40.2357, 2.57338],
        errors=[288.248, 202.527, 134.025, 27.1269, 12.3450, 1.13894, 0.213825],
    ),
    # ── DEM 99 ────────────────────────────────────────────────────────────────
    dict(
        dem_id=99,
        nfilters=4,
        filters=["Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[35197.6, 8244.95, 3381.95, 338.260],
        errors=[1063.40, 261.701, 99.4482, 9.33411],
    ),
    dict(
        dem_id=99,
        nfilters=5,
        filters=["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[59826.3, 35197.6, 8244.95, 3381.95, 338.260],
        errors=[1620.93, 1063.40, 261.701, 99.4482, 9.33411],
    ),
    dict(
        dem_id=99,
        nfilters=6,
        filters=["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[87662.9, 59826.3, 35197.6, 8244.95, 3381.95, 338.260],
        errors=[2403.84, 1620.93, 1063.40, 261.701, 99.4482, 9.33411],
    ),
    dict(
        dem_id=99,
        nfilters=7,
        filters=[
            "C-poly",
            "Ti-poly",
            "Be-thin",
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
        ],
        intensities=[87662.9, 59826.3, 35197.6, 8244.95, 3381.95, 338.260, 22.4599],
        errors=[2403.84, 1620.93, 1063.40, 261.701, 99.4482, 9.33411, 0.718465],
    ),
    # ── DEM 112 ───────────────────────────────────────────────────────────────
    dict(
        dem_id=112,
        nfilters=4,
        filters=["Be-med", "Al-med", "Al-thick", "Be-thick"],
        intensities=[215903.0, 93715.4, 8685.49, 330.326],
        errors=[6773.48, 2878.83, 268.391, 9.22694],
    ),
    dict(
        dem_id=112,
        nfilters=5,
        filters=["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"],
        intensities=[215903.0, 93715.4, 8685.49, 330.326, 8099.61],
        errors=[6773.48, 2878.83, 268.391, 9.22694, 259.455],
    ),
    dict(
        dem_id=112,
        nfilters=6,
        filters=[
            "Al-thick",
            "Be-thick",
            "Al-poly/Al-mesh",
            "Al-poly/Al-thick",
            "Al-poly/Be-thick",
            "C-poly/Ti-poly",
        ],
        intensities=[8685.49, 330.326, 3295230.0, 8099.61, 293.908, 1253340.0],
        errors=[268.391, 9.22694, 58410.6, 259.455, 8.71760, 36413.1],
    ),
    dict(
        dem_id=112,
        nfilters=7,
        filters=[
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
            "Al-poly/Al-thick",
            "Al-poly/Be-thick",
            "C-poly/Ti-poly",
        ],
        intensities=[215903.0, 93715.4, 8685.49, 330.326, 8099.61, 293.908, 1253340.0],
        errors=[6773.48, 2878.83, 268.391, 9.22694, 259.455, 8.71760, 36413.1],
    ),
    # ── DEM 159 ───────────────────────────────────────────────────────────────
    dict(
        dem_id=159,
        nfilters=4,
        filters=["Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[5438.57, 1621.40, 654.170, 72.3524],
        errors=[149.700, 48.7544, 18.3656, 2.16314],
    ),
    dict(
        dem_id=159,
        nfilters=5,
        filters=["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[8086.46, 5438.57, 1621.40, 654.170, 72.3524],
        errors=[216.947, 149.700, 48.7544, 18.3656, 2.16314],
    ),
    dict(
        dem_id=159,
        nfilters=6,
        filters=["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities=[11054.7, 8086.46, 5438.57, 1621.40, 654.170, 72.3524],
        errors=[301.769, 216.947, 149.700, 48.7544, 18.3656, 2.16314],
    ),
    dict(
        dem_id=159,
        nfilters=7,
        filters=[
            "C-poly",
            "Ti-poly",
            "Be-thin",
            "Be-med",
            "Al-med",
            "Al-thick",
            "Be-thick",
        ],
        intensities=[11054.7, 8086.46, 5438.57, 1621.40, 654.170, 72.3524, 6.14447],
        errors=[301.769, 216.947, 149.700, 48.7544, 18.3656, 2.16314, 0.321695],
    ),
]


def load_true_dem(dem_id):
    """Load true DEM — returns logT (26,), dem (26,)."""
    txt = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    data = np.loadtxt(txt)
    return data[:, 0], data[:, 1]


def load_idl_sav(dem_id, nfilters):
    """Load IDL base DEM from .sav with errors — returns logT, dem (26,), chisq."""
    sav = (
        IDL_SAV_DIR
        / f"idl_dem_idx{dem_id}_base_wErrors_10Oct2008_{nfilters}filters.sav"
    )
    data = readsav(str(sav), python_dict=True)
    logT = np.array(data["logt_out"]).ravel()
    dem = np.array(data["dem_out"]).ravel()
    chisq = float(data["chisq"])
    return logT, dem, chisq


def run_xrtpy(filters, intensities, errors):
    """Run XRTpy base DEM with errors — returns logT, dem, chisq."""
    responses = generate_temperature_responses(filters, OBSERVATION_DATE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver = XRTDEMIterative(
            observed_channel=filters,
            observed_intensities=intensities,
            temperature_responses=responses,
            intensity_uncertainties=errors,
            monte_carlo_runs=0,
        )
        solver.solve()
    return solver.logT, solver.dem, solver.chisq


def run_and_plot(run):
    dem_id = run["dem_id"]
    nfilters = run["nfilters"]
    filters = run["filters"]
    intensities = run["intensities"]
    errors = run["errors"]

    print(f"  DEM {dem_id} — {nfilters} filters")

    # ── Load all three ────────────────────────────────────────────────────────
    logT_true, dem_true = load_true_dem(dem_id)
    logT_idl, dem_idl, chisq_idl = load_idl_sav(dem_id, nfilters)
    logT_py, dem_py, chisq_py = run_xrtpy(filters, intensities, errors)

    log_true = np.log10(np.clip(dem_true, 1e-40, None))
    log_idl = np.log10(np.clip(dem_idl, 1e-40, None))
    log_py = np.log10(np.clip(dem_py, 1e-40, None))

    # ── Filter label (intensity ± error) ─────────────────────────────────────
    filter_label = ",  ".join(
        f"{f} [{i:.4g}±{e:.3g}]" for f, i, e in zip(filters, intensities, errors)
    )

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))

    ax.step(
        logT_true, log_true, where="mid", color="black", linewidth=2.5, label="True DEM"
    )

    ax.step(
        logT_idl,
        log_idl,
        where="mid",
        color="orange",
        linewidth=2.0,
        label=f"IDL  |  χ² = {chisq_idl:.4f}",
    )

    ax.step(
        logT_py,
        log_py,
        where="mid",
        color="#1E90FF",
        linewidth=2.0,
        linestyle="--",
        label=f"XRTpy  |  χ² = {chisq_py:.4f}",
    )

    ax.set_title(
        f"Hinode/XRT — Base DEM Comparison w/ Errors — DEM Index {dem_id}",
        fontsize=13,
        pad=12,
    )
    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

    fig.text(
        0.5, 0.91, filter_label, ha="center", va="top", fontsize=7.5, color="black"
    )

    ax.legend(fontsize=11)
    ax.grid(visible=True, alpha=0.3)
    ax.set_xlim(5.5, 8.0)
    ax.set_ylim(14, 26)

    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        OUT_DIR / f"dem_idx{dem_id}_base_wErrors_IDL_XRTpy_True_"
        f"{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_"
        f"{nfilters}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"    Saved: {outfile.name}")


# ── Run all 24 ────────────────────────────────────────────────────────────────
print("Running all 24 base DEM comparisons with errors...")
for run in RUNS:
    run_and_plot(run)

print("\nAll done — 24 plots saved to xrtpy_idl_fake_data_base_dem_comparison/")

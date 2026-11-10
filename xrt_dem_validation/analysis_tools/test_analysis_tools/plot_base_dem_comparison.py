# plot_base_dem_comparison.py
#
# For each of the 14 runs plots THREE curves on one figure:
#   - Orange solid  : IDL base DEM  (from MC0 .sav)
#   - Blue dashed   : XRTpy base DEM (solved fresh)
#   - Black solid   : True DEM (from DEM_XX.txt)

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR     = Path("/Users/jvelasq/Projects/xrtpy/xrtpy/xrt_dem_iterative")
MC0_DIR      = BASE_DIR / "Fake_IDL_data"
TRUE_DEM_DIR = BASE_DIR / "Fake_IDL_data" / "Fake_Based_IDL_DEMs"
OUT_DIR      = BASE_DIR  # plots saved here

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL       = "10-Oct-2008 00:00:00"

# ── All 14 runs ───────────────────────────────────────────────────────────────
RUNS = [
    dict(dem_id=14,  nfilters=7,
         filters     = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [4949.20,   5291.10,   2329.27,   1112.14,  216.758,  94.3998,  8.71282]),
    dict(dem_id=15,  nfilters=7,
         filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
         intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174,   1.48818]),
    dict(dem_id=15,  nfilters=6,
         filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174]),
    dict(dem_id=65,  nfilters=5,
         filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [7020.42,   4331.36,  1007.78,   416.894,  40.2357]),
    dict(dem_id=65,  nfilters=6,
         filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [10060.6,   7020.42,  4331.36,   1007.78,  416.894,  40.2357]),
    dict(dem_id=99,  nfilters=5,
         filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [59826.3,   35197.6,  8244.95,   3381.95,  338.260]),
    dict(dem_id=99,  nfilters=6,
         filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [87662.9,   59826.3,  35197.6,   8244.95,  3381.95,  338.260]),
    dict(dem_id=99,  nfilters=7,
         filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
         intensities = [87662.9,   59826.3,  35197.6,   8244.95,  3381.95,  338.260,   22.4599]),
    dict(dem_id=112, nfilters=5,
         filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"],
         intensities = [215903.,   93715.4,  8685.49,   330.326,   8099.61]),
    dict(dem_id=112, nfilters=6,
         filters     = ["Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
         intensities = [8685.49,    330.326,   3295230.,            8099.61,            293.908,            1253340.]),
    dict(dem_id=112, nfilters=7,
         filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
         intensities = [215903.,   93715.4,  8685.49,   330.326,   8099.61,            293.908,            1253340.]),
    dict(dem_id=159, nfilters=5,
         filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [8086.46,   5438.57,  1621.40,   654.170,  72.3524]),
    dict(dem_id=159, nfilters=6,
         filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
         intensities = [11054.7,   8086.46,  5438.57,   1621.40,  654.170,  72.3524]),
    dict(dem_id=159, nfilters=7,
         filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
         intensities = [11054.7,   8086.46,  5438.57,   1621.40,  654.170,  72.3524,   6.14447]),
]


def load_true_dem(dem_id):
    """Load true DEM from DEM_XX.txt — returns logT (26,), dem (26,)."""
    txt = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    data = np.loadtxt(txt)
    return data[:, 0], data[:, 1]


def load_idl_base(dem_id, nfilters):
    """Load IDL base DEM from MC0 .sav — returns logT, dem (26,), chisq."""
    sav = MC0_DIR / f"idl_dem_idx{dem_id}_10Oct2008_{nfilters}filters_MC0.sav"
    data = readsav(str(sav), python_dict=True)
    logT  = np.array(data['logt_out']).ravel()
    dem   = np.array(data['dem_out']).ravel()
    chisq = float(data['chisq'])
    return logT, dem, chisq


def run_xrtpy_base(filters, intensities):
    """Run XRTpy MC=0 — returns logT, dem, chisq."""
    responses = generate_temperature_responses(filters, OBSERVATION_DATE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver = XRTDEMIterative(
            observed_channel=filters,
            observed_intensities=intensities,
            temperature_responses=responses,
            monte_carlo_runs=0,
        )
        solver.solve()
    return solver.logT, solver.dem, solver.chisq


def run_and_plot(run):
    dem_id   = run['dem_id']
    nfilters = run['nfilters']
    filters  = run['filters']
    intensities = run['intensities']

    print(f"  DEM {dem_id} — {nfilters} filters")

    # ── Load all three DEMs ───────────────────────────────────────────────────
    logT_true, dem_true   = load_true_dem(dem_id)
    logT_idl,  dem_idl,  chisq_idl = load_idl_base(dem_id, nfilters)
    logT_py,   dem_py,   chisq_py  = run_xrtpy_base(filters, intensities)

    # ── log10 ─────────────────────────────────────────────────────────────────
    log_true = np.log10(np.clip(dem_true, 1e-40, None))
    log_idl  = np.log10(np.clip(dem_idl,  1e-40, None))
    log_py   = np.log10(np.clip(dem_py,   1e-40, None))

    # ── Filter label ──────────────────────────────────────────────────────────
    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    # ── Plot ──────────────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(10, 6))

    # True DEM — black solid, drawn first (background reference)
    ax.step(logT_true, log_true, where="mid",
            color="black", linewidth=2.5,
            label="True DEM")

    # IDL base — orange solid
    ax.step(logT_idl, log_idl, where="mid",
            color="orange", linewidth=2.0,
            label=f"IDL  |  χ² = {chisq_idl:.4f}")

    # XRTpy base — Dodger Blue dashed
    ax.step(logT_py, log_py, where="mid",
            color="#1E90FF", linewidth=2.0, linestyle="--",
            label=f"XRTpy  |  χ² = {chisq_py:.4f}")

    ax.set_title(
        f"Hinode/XRT — Base DEM Comparison — DEM Index {dem_id}",
        fontsize=13, pad=12,
    )
    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)
    fig.text(0.5, 0.91, filter_label, ha="center", va="top",
             fontsize=8.5, color="black")

    ax.legend(fontsize=11)
    ax.grid(visible=True, alpha=0.3)
    ax.set_xlim(5.5, 8.0)
    #ax.set_xlim(logT_py.min(), logT_py.max())
    #ax.set_ylim(5.5, 8.0)
    ax.set_ylim(16, 26)

     # y-limits from all three curves

     #     all_log = np.concatenate([log_true, log_idl, log_py])
     #     finite  = np.isfinite(all_log)
     #     if np.any(finite):
     #         y   = all_log[finite]
     #         pad = 0.25 * (y.max() - y.min() + 1e-6)
     #         ax.set_ylim(y.min() - pad, y.max() + pad)


    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        OUT_DIR /
        f"dem_idx{dem_id}_base_IDL_XRTpy_True_"
        f"{DATE_LABEL.replace(' ','_').replace(':','').replace('-','')}_"
        f"{nfilters}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"    Saved: {outfile.name}")


# ── Run all ───────────────────────────────────────────────────────────────────
print("Running all 14 base DEM comparisons...")
for run in RUNS:
    run_and_plot(run)

print("\nAll done — 14 plots saved.")

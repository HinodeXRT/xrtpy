# plot_mc100_comparison.py
#
# For each of the 14 runs produces TWO plots:
#   1) IDL MC100 only  (base black + 100 grey MC curves)
#   2) IDL MC100 + XRTpy MC100 overlaid
#      - IDL base : orange solid
#      - IDL MC   : orange, thin, alpha
#      - XRTpy base: Dodger Blue dashed
#      - XRTpy MC  : Dodger Blue, thin, alpha

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ── Paths ─────────────────────────────────────────────────────────────────────
BASE_DIR    = Path("/Users/jvelasq/Projects/xrtpy/xrtpy/xrt_dem_iterative")
MC0_DIR     = BASE_DIR / "Fake_IDL_data"
MC100_DIR   = BASE_DIR / "Fake_IDL_data" / "Fake_IDL_MC_data"

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL       = "10-Oct-2008 00:00:00"
MC_RUNS          = 100

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


def load_idl_mc(dem_id, nfilters):
    """Load IDL MC100 .sav — returns logT, dem_base (26,), dem_mc (100,26), chisq_base."""
    sav = MC100_DIR / f"idl_dem_idx{dem_id}_10Oct2008_{nfilters}filters_MC100.sav"
    data = readsav(str(sav), python_dict=True)

    logT = np.array(data['logt_out']).ravel()
    dem_raw = np.array(data['dem_out'])

    # orient to (nT, nRuns+1) then transpose to (nRuns+1, nT)
    nT = len(logT)
    if dem_raw.ndim == 1:
        # MC0 accidentally saved — treat as base only
        dem_all = dem_raw.reshape(1, nT)
    elif dem_raw.shape[0] == nT:
        dem_all = dem_raw.T          # was (nT, nRuns+1)
    else:
        dem_all = dem_raw            # already (nRuns+1, nT)

    dem_base = dem_all[0]            # (nT,)
    dem_mc   = dem_all[1:]           # (100, nT)
    chisq    = float(data['chisq']) if np.array(data['chisq']).ndim == 0 else float(data['chisq'][0])

    return logT, dem_base, dem_mc, chisq


def run_xrtpy_mc(filters, intensities):
    """Run XRTpy with MC=100, return logT, dem_base, dem_mc (100,nT), chisq."""
    responses = generate_temperature_responses(filters, OBSERVATION_DATE)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        solver = XRTDEMIterative(
            observed_channel=filters,
            observed_intensities=intensities,
            temperature_responses=responses,
            monte_carlo_runs=MC_RUNS,
        )
        solver.solve()

    logT     = solver.logT
    dem_base = solver.dem                 # (nT,)
    dem_mc   = solver.mc_dem[1:, :]      # (100, nT)  rows 1..100
    chisq    = solver.chisq

    return logT, dem_base, dem_mc, chisq


def _log10(arr):
    return np.log10(np.clip(arr, 1e-40, None))


def _ylim(arrays):
    """Compute y-limits from a list of log10 arrays."""
    all_vals = np.concatenate([a.ravel() for a in arrays])
    finite = all_vals[np.isfinite(all_vals)]
    if len(finite) == 0:
        return -5, 30
    pad = 0.15 * (finite.max() - finite.min() + 1e-6)
    return finite.min() - pad, finite.max() + pad


def plot_idl_mc_only(dem_id, nfilters, filters, intensities,
                     logT_idl, dem_base_idl, dem_mc_idl, chisq_idl):
    """Plot 1: IDL base DEM + 100 MC curves only."""

    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    log_base = _log10(dem_base_idl)
    log_mc   = _log10(dem_mc_idl)          # (100, nT)

    fig, ax = plt.subplots(figsize=(10, 6))

    # MC curves first
    for i in range(log_mc.shape[0]):
        ax.step(logT_idl, log_mc[i], where="mid",
                color="black", linewidth=0.8, alpha=0.18)

    # Base on top
    ax.step(logT_idl, log_base, where="mid",
            color="black", linewidth=2.5,
            label=f"{DATE_LABEL}  |  χ² = {chisq_idl:.4f}")

    ax.set_title(
        f"Hinode/XRT/IDL — Base DEM + {MC_RUNS} MC — DEM Index {dem_id}",
        fontsize=13, pad=12,
    )
    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)
    fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

    ax.legend(fontsize=11)
    ax.grid(visible=True, alpha=0.3)
    ax.set_xlim(logT_idl.min(), logT_idl.max())
    ax.set_ylim(*_ylim([log_base, log_mc.ravel()]))

    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        BASE_DIR /
        f"dem_idx{dem_id}_IDL_MC{MC_RUNS}_"
        f"{DATE_LABEL.replace(' ','_').replace(':','').replace('-','')}_"
        f"{nfilters}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outfile.name}")


def plot_idl_vs_xrtpy_mc(dem_id, nfilters, filters, intensities,
                         logT_idl, dem_base_idl, dem_mc_idl, chisq_idl,
                         logT_py,  dem_base_py,  dem_mc_py,  chisq_py):
    """Plot 2: IDL MC100 + XRTpy MC100 overlaid."""

    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    log_base_idl = _log10(dem_base_idl)
    log_mc_idl   = _log10(dem_mc_idl)
    log_base_py  = _log10(dem_base_py)
    log_mc_py    = _log10(dem_mc_py)

    fig, ax = plt.subplots(figsize=(10, 6))

    # IDL MC (orange, thin)
    for i in range(log_mc_idl.shape[0]):
        ax.step(logT_idl, log_mc_idl[i], where="mid",
                color="orange", linewidth=0.8, alpha=0.18)

    # XRTpy MC (blue, thin)
    for i in range(log_mc_py.shape[0]):
        ax.step(logT_py, log_mc_py[i], where="mid",
                color="#1E90FF", linewidth=0.8, alpha=0.18)

    # IDL base (orange, thick solid)
    ax.step(logT_idl, log_base_idl, where="mid",
            color="orange", linewidth=2.5,
            label=f"IDL  |  χ² = {chisq_idl:.4f}")

    # XRTpy base (blue, thick dashed)
    ax.step(logT_py, log_base_py, where="mid",
            color="#1E90FF", linewidth=2.5, linestyle="--",
            label=f"XRTpy  |  χ² = {chisq_py:.4f}")

    ax.set_title(
        f"Hinode/XRT — IDL vs XRTpy — Base + {MC_RUNS} MC — DEM Index {dem_id}",
        fontsize=13, pad=12,
    )
    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)
    fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

    ax.legend(fontsize=11)
    ax.grid(visible=True, alpha=0.3)
    ax.set_xlim(logT_py.min(), logT_py.max())
    ax.set_ylim(*_ylim([log_base_idl, log_mc_idl.ravel(),
                        log_base_py,  log_mc_py.ravel()]))

    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        BASE_DIR /
        f"dem_idx{dem_id}_IDLvsXRTpy_MC{MC_RUNS}_"
        f"{DATE_LABEL.replace(' ','_').replace(':','').replace('-','')}_"
        f"{nfilters}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"  Saved: {outfile.name}")


# ── Main loop ─────────────────────────────────────────────────────────────────
for run in RUNS:
    dem_id   = run['dem_id']
    nfilters = run['nfilters']
    filters  = run['filters']
    intensities = run['intensities']

    print(f"\nDEM {dem_id} — {nfilters} filters")

    # Load IDL MC100
    logT_idl, dem_base_idl, dem_mc_idl, chisq_idl = load_idl_mc(dem_id, nfilters)

    # Run XRTpy MC100
    print(f"  Running XRTpy MC{MC_RUNS}...")
    logT_py, dem_base_py, dem_mc_py, chisq_py = run_xrtpy_mc(filters, intensities)

    # Plot 1: IDL only
    plot_idl_mc_only(dem_id, nfilters, filters, intensities,
                     logT_idl, dem_base_idl, dem_mc_idl, chisq_idl)

    # Plot 2: IDL + XRTpy overlaid
    plot_idl_vs_xrtpy_mc(dem_id, nfilters, filters, intensities,
                         logT_idl, dem_base_idl, dem_mc_idl, chisq_idl,
                         logT_py,  dem_base_py,  dem_mc_py,  chisq_py)

print("\nAll done — 28 plots saved.")

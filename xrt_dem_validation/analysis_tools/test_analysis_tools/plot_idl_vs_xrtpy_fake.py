# plot_all_idl_vs_xrtpy.py

import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL       = "10-Oct-2008 00:00:00"
SAV_DIR          = Path("Fake_IDL_data")

# ── All 14 runs ───────────────────────────────────────────────────────────────
RUNS = [
    dict(
        dem_id   = 14,
        nfilters = 7,
        filters     = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [4949.20,   5291.10,   2329.27,   1112.14,  216.758,  94.3998,  8.71282],
    ),
    dict(
        dem_id   = 15,
        nfilters = 7,
        filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
        intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174,   1.48818],
    ),
    dict(
        dem_id   = 15,
        nfilters = 6,
        filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174],
    ),
    dict(
        dem_id   = 65,
        nfilters = 5,
        filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [7020.42,   4331.36,  1007.78,   416.894,  40.2357],
    ),
    dict(
        dem_id   = 65,
        nfilters = 6,
        filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [10060.6,   7020.42,  4331.36,   1007.78,  416.894,  40.2357],
    ),
    dict(
        dem_id   = 99,
        nfilters = 5,
        filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [59826.3,   35197.6,  8244.95,   3381.95,  338.260],
    ),
    dict(
        dem_id   = 99,
        nfilters = 6,
        filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [87662.9,   59826.3,  35197.6,   8244.95,  3381.95,  338.260],
    ),
    dict(
        dem_id   = 99,
        nfilters = 7,
        filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
        intensities = [87662.9,   59826.3,  35197.6,   8244.95,  3381.95,  338.260,   22.4599],
    ),
    dict(
        dem_id   = 112,
        nfilters = 5,
        filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"],
        intensities = [215903.,   93715.4,  8685.49,   330.326,   8099.61],
    ),
    dict(
        dem_id   = 112,
        nfilters = 6,
        filters     = ["Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
        intensities = [8685.49,    330.326,   3295230.,            8099.61,            293.908,            1253340.],
    ),
    dict(
        dem_id   = 112,
        nfilters = 7,
        filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
        intensities = [215903.,   93715.4,  8685.49,   330.326,   8099.61,            293.908,            1253340.],
    ),
    dict(
        dem_id   = 159,
        nfilters = 5,
        filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [8086.46,   5438.57,  1621.40,   654.170,  72.3524],
    ),
    dict(
        dem_id   = 159,
        nfilters = 6,
        filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
        intensities = [11054.7,   8086.46,  5438.57,   1621.40,  654.170,  72.3524],
    ),
    dict(
        dem_id   = 159,
        nfilters = 7,
        filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
        intensities = [11054.7,   8086.46,  5438.57,   1621.40,  654.170,  72.3524,   6.14447],
    ),
]


def run_and_plot(run):
    dem_id   = run['dem_id']
    nfilters = run['nfilters']
    filters  = run['filters']
    intensities = run['intensities']

    sav_path = SAV_DIR / f"idl_dem_idx{dem_id}_10Oct2008_{nfilters}filters_MC0.sav"

    # ── Load IDL .sav ─────────────────────────────────────────────────────────
    data     = readsav(str(sav_path), python_dict=True)
    logT_idl = np.array(data['logt_out']).ravel()
    dem_idl  = np.array(data['dem_out']).ravel()
    chisq_idl = float(data['chisq'])

    # ── Run XRTpy ────────────────────────────────────────────────────────────
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

    logT_py   = solver.logT
    dem_py    = solver.dem
    chisq_py  = solver.chisq

    # ── Sanity check ─────────────────────────────────────────────────────────
    assert np.allclose(logT_idl, logT_py, atol=1e-3), \
        f"logT grids differ for DEM {dem_id} {nfilters}f"

    # ── Plot ─────────────────────────────────────────────────────────────────
    log_dem_idl = np.log10(np.clip(dem_idl, 1e-40, None))
    log_dem_py  = np.log10(np.clip(dem_py,  1e-40, None))

    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.step(logT_idl, log_dem_idl, where="mid",
            color="orange", linewidth=2.5,
            label=f"IDL  |  χ² = {chisq_idl:.4f}")

    ax.step(logT_py, log_dem_py, where="mid",
            color="#1E90FF", linewidth=2.5, linestyle="--",
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
    ax.set_xlim(logT_py.min(), logT_py.max())

    all_log = np.concatenate([log_dem_idl, log_dem_py])
    finite  = np.isfinite(all_log)
    if np.any(finite):
        y   = all_log[finite]
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        ax.set_ylim(y.min() - pad, y.max() + pad)

    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        f"dem_idx{dem_id}_IDLvsXRTpy_"
        f"{DATE_LABEL.replace(' ','_').replace(':','').replace('-','')}_"
        f"{nfilters}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.close(fig)
    print(f"Saved: {outfile}")


# ── Run all ───────────────────────────────────────────────────────────────────
for run in RUNS:
    run_and_plot(run)

print("\nAll 14 comparison plots done.")

# import matplotlib.pyplot as plt
# import numpy as np
# from pathlib import Path
# from scipy.io import readsav
# from xrtpy.response.tools import generate_temperature_responses
# from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # ── CONFIG — change these to match whichever .sav you are comparing ──────────
# # SAV_PATH  = Path("Fake_IDL_data/idl_dem_idx159_10Oct2008_7filters_MC0.sav")
# # DEM_ID    = 159
# # DATE_LABEL = "10-Oct-2008 00:00:00"
# # OBSERVATION_DATE = "2008-10-10T00:00:00"

# # # Must match exactly what is in the .sav (same order)
# # filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"]
# # intensities = [11054.7,   8086.46,   5438.57,  1621.40,  654.170,  72.3524,   6.14447]

# # # ── CONFIG ───────────────────────────────────────────────────────────────────
# # SAV_PATH  = Path("Fake_IDL_data/idl_dem_idx14_10Oct2008_7filters_MC0.sav")
# # DEM_ID    = 14
# # DATE_LABEL = "10-Oct-2008 00:00:00"
# # OBSERVATION_DATE = "2008-10-10T00:00:00"

# # filters     = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# # intensities = [4949.20,   5291.10,   2329.27,   1112.14,  216.758,  94.3998,  8.71282]

# # ── Load IDL .sav ─────────────────────────────────────────────────────────────
# data    = readsav(str(SAV_PATH), python_dict=True)
# logT_idl = np.array(data['logt_out']).ravel()
# dem_idl  = np.array(data['dem_out']).ravel()
# chisq_idl = float(data['chisq'])

# # ── Run XRTpy DEM ─────────────────────────────────────────────────────────────
# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# solver.solve()

# logT_py  = solver.logT
# dem_py   = solver.dem
# chisq_py = solver.chisq

# # ── Sanity check grids match ──────────────────────────────────────────────────
# assert np.allclose(logT_idl, logT_py, atol=1e-3), \
#     f"logT grids differ! IDL: {logT_idl[:3]}  XRTpy: {logT_py[:3]}"

# # ── Plot ──────────────────────────────────────────────────────────────────────
# log_dem_idl = np.log10(np.clip(dem_idl, 1e-40, None))
# log_dem_py  = np.log10(np.clip(dem_py,  1e-40, None))

# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT_idl, log_dem_idl, where="mid",
#         color="orange", linewidth=2.5,
#         label=f"IDL  |  χ² = {chisq_idl:.4f}")

# ax.step(logT_py, log_dem_py, where="mid",
#         color="#1E90FF", linewidth=2.5, linestyle="--",
#         label=f"XRTpy  |  χ² = {chisq_py:.4f}")

# ax.set_title(
#     f"Hinode/XRT — Base DEM Comparison — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT_py.min(), logT_py.max())

# # y-limits from both curves combined
# all_log = np.concatenate([log_dem_idl, log_dem_py])
# finite  = np.isfinite(all_log)
# if np.any(finite):
#     y   = all_log[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_IDLvsXRTpy_{DATE_LABEL.replace(' ','_').replace(':','').replace('-','')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")

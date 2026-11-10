import matplotlib.pyplot as plt
import numpy as np

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL       = "10-Oct-2008 00:00:00"
MC_RUNS          = 100


# ── Shared plot + save helper ─────────────────────────────────────────────────
def run_and_plot_mc(DEM_ID, filters, intensities):

    responses = generate_temperature_responses(filters, OBSERVATION_DATE)

    dem_solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=MC_RUNS,
    )
    dem_solver.solve()

    dem    = dem_solver.dem
    logT   = dem_solver.logT
    chisq  = dem_solver.chisq
    mc_dem = dem_solver.mc_dem  # shape (MC_RUNS+1, nT) — row 0 = base

    dem_safe = np.clip(dem, 1e-40, None)
    log_dem  = np.log10(dem_safe)

    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    fig, ax = plt.subplots(figsize=(10, 6))

    # MC curves — black, thin, transparent — drawn first so base sits on top
    for i in range(1, mc_dem.shape[0]):
        mc_safe = np.clip(mc_dem[i], 1e-40, None)
        ax.step(logT, np.log10(mc_safe), where="mid",
                color="black", linewidth=0.8, alpha=0.18)

    # Base DEM — Dodger Blue on top
    ax.step(logT, log_dem, where="mid",
            color="#1E90FF", linewidth=2.5,
            label=f"{DATE_LABEL}  |  χ² = {chisq:.4f}")

    ax.set_title(
        f"Hinode/XRT/XRTpy — Base DEM + {MC_RUNS} MC — DEM Index {DEM_ID}",
        fontsize=13, pad=12,
    )
    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

    fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

    ax.legend(fontsize=11)
    ax.grid(visible=True, alpha=0.3)
    ax.set_xlim(logT.min(), logT.max())

    finite = np.isfinite(log_dem)
    if np.any(finite):
        y   = log_dem[finite]
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        ax.set_ylim(y.min() - pad, y.max() + pad)

    plt.tight_layout(rect=[0, 0, 1, 0.91])

    outfile = (
        f"dem_idx{DEM_ID}_"
        f"{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_"
        f"{len(filters)}filters_MC{MC_RUNS}.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.show()
    print(f"Saved: {outfile}")


##################################################################################################################################
##################################      INDEX 14                              ####################################################
##################################################################################################################################

run_and_plot_mc(
    DEM_ID      = 14,
    filters     = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [4949.20,   5291.10,   2329.27,   1112.14,  216.758,  94.3998,  8.71282],
)

##################################################################################################################################
##################################      INDEX 15                              ####################################################
##################################################################################################################################

# 7 filters
run_and_plot_mc(
    DEM_ID      = 15,
    filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174,   1.48818],
)

# 6 filters
run_and_plot_mc(
    DEM_ID      = 15,
    filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174],
)

##################################################################################################################################
##################################      INDEX 65                              ####################################################
##################################################################################################################################

# 5 filters
run_and_plot_mc(
    DEM_ID      = 65,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [7020.42,   4331.36,   1007.78,  416.894,  40.2357],
)

# 6 filters
run_and_plot_mc(
    DEM_ID      = 65,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [10060.6,  7020.42,   4331.36,   1007.78,  416.894,  40.2357],
)

##################################################################################################################################
##################################      INDEX 99                              ####################################################
##################################################################################################################################

# 5 filters
run_and_plot_mc(
    DEM_ID      = 99,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [59826.3,   35197.6,   8244.95,  3381.95,  338.260],
)

# 6 filters
run_and_plot_mc(
    DEM_ID      = 99,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [87662.9,  59826.3,   35197.6,   8244.95,  3381.95,  338.260],
)

# 7 filters
run_and_plot_mc(
    DEM_ID      = 99,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [87662.9,  59826.3,   35197.6,   8244.95,  3381.95,  338.260,   22.4599],
)

##################################################################################################################################
##################################      INDEX 112                             ####################################################
##################################################################################################################################

# 5 filters
run_and_plot_mc(
    DEM_ID      = 112,
    filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"],
    intensities = [215903.,  93715.4,  8685.49,   330.326,   8099.61],
)

# 6 filters
run_and_plot_mc(
    DEM_ID      = 112,
    filters     = ["Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
    intensities = [8685.49,   330.326,   3295230.,            8099.61,            293.908,            1253340.],
)

# 7 filters
run_and_plot_mc(
    DEM_ID      = 112,
    filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
    intensities = [215903.,  93715.4,  8685.49,   330.326,   8099.61,            293.908,            1253340.],
)

##################################################################################################################################
##################################      INDEX 159                             ####################################################
##################################################################################################################################

# 5 filters
run_and_plot_mc(
    DEM_ID      = 159,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [8086.46,   5438.57,   1621.40,  654.170,  72.3524],
)

# 6 filters
run_and_plot_mc(
    DEM_ID      = 159,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [11054.7,  8086.46,   5438.57,   1621.40,  654.170,  72.3524],
)

# 7 filters
run_and_plot_mc(
    DEM_ID      = 159,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [11054.7,  8086.46,   5438.57,   1621.40,  654.170,  72.3524,   6.14447],
)

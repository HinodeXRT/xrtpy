import matplotlib.pyplot as plt
import numpy as np

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

OBSERVATION_DATE = "2008-10-10T00:00:00"
DATE_LABEL       = "10-Oct-2008 00:00:00"


# ── Shared plot + save helper ─────────────────────────────────────────────────
def run_and_plot(DEM_ID, filters, intensities):

    responses = generate_temperature_responses(filters, OBSERVATION_DATE)

    dem_solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=0,
    )
    dem_solver.solve()

    dem   = dem_solver.dem
    logT  = dem_solver.logT
    chisq = dem_solver.chisq

    dem_safe     = np.clip(dem, 1e-40, None)
    log_dem      = np.log10(dem_safe)
    filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0,
            label=f"{DATE_LABEL}  |  χ² = {chisq:.4f}")

    ax.set_title(
        f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
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
        f"{len(filters)}filters.png"
    )
    fig.savefig(outfile, dpi=150)
    plt.show()
    print(f"Saved: {outfile}")


##################################################################################################################################
##################################      INDEX 14                              ####################################################
##################################################################################################################################

run_and_plot(
    DEM_ID      = 14,
    filters     = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [4949.20,   5291.10,   2329.27,   1112.14,  216.758,  94.3998,  8.71282],
)

##################################################################################################################################
##################################      INDEX 15                              ####################################################
##################################################################################################################################

# 7 filters
run_and_plot(
    DEM_ID      = 15,
    filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174,   1.48818],
)

# 6 filters
run_and_plot(
    DEM_ID      = 15,
    filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174],
)

##################################################################################################################################
##################################      INDEX 65                              ####################################################
##################################################################################################################################

# 5 filters
run_and_plot(
    DEM_ID      = 65,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [7020.42,   4331.36,   1007.78,  416.894,  40.2357],
)

# 6 filters
run_and_plot(
    DEM_ID      = 65,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [10060.6,   7020.42,   4331.36,  1007.78,  416.894,  40.2357],
)

##################################################################################################################################
##################################      INDEX 99                              ####################################################
##################################################################################################################################

# 5 filters
run_and_plot(
    DEM_ID      = 99,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [59826.3,   35197.6,   8244.95,  3381.95,  338.260],
)

# 6 filters
run_and_plot(
    DEM_ID      = 99,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [87662.9,   59826.3,   35197.6,  8244.95,  3381.95,  338.260],
)

# 7 filters
run_and_plot(
    DEM_ID      = 99,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [87662.9,   59826.3,   35197.6,  8244.95,  3381.95,  338.260,   22.4599],
)

##################################################################################################################################
##################################      INDEX 112                             ####################################################
##################################################################################################################################

# 5 filters
run_and_plot(
    DEM_ID      = 112,
    filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"],
    intensities = [215903.,  93715.4,  8685.49,    330.326,   8099.61],
)

# 6 filters
run_and_plot(
    DEM_ID      = 112,
    filters     = ["Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
    intensities = [8685.49,    330.326,   3295230.,            8099.61,            293.908,            1253340.],
)

# 7 filters
run_and_plot(
    DEM_ID      = 112,
    filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"],
    intensities = [215903.,  93715.4,  8685.49,    330.326,   8099.61,            293.908,            1253340.],
)

##################################################################################################################################
##################################      INDEX 159                             ####################################################
##################################################################################################################################

# 5 filters
run_and_plot(
    DEM_ID      = 159,
    filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [8086.46,   5438.57,   1621.40,  654.170,  72.3524],
)

# 6 filters
run_and_plot(
    DEM_ID      = 159,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"],
    intensities = [11054.7,   8086.46,   5438.57,  1621.40,  654.170,  72.3524],
)

# 7 filters
run_and_plot(
    DEM_ID      = 159,
    filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"],
    intensities = [11054.7,   8086.46,   5438.57,  1621.40,  654.170,  72.3524,   6.14447],
)

# # from xrtpy.response.tools import generate_temperature_responses
# # from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # # DEM 14 — 2008-10-10
# # filters = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# # intensities = [4949.20, 5291.10, 2329.27, 1112.14, 216.758, 94.3998, 8.71282]
# # observation_date = "2008-10-10T00:00:00"

# # responses = generate_temperature_responses(filters, observation_date)

# # dem_solver = XRTDEMIterative(
# #     observed_channel=filters,
# #     observed_intensities=intensities,
# #     temperature_responses=responses,
# #     monte_carlo_runs=0,
# # )

# # dem_solver.solve()
# # dem_solver.plot_dem()

# import matplotlib.pyplot as plt
# import numpy as np
# from xrtpy.response.tools import generate_temperature_responses
# from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # ── DEM 14 metadata ──────────────────────────────────────────────────────────
# DEM_ID          = 14
# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL      = "10-Oct-2008 00:00:00"

# filters = ["Al-mesh", "Al-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [4949.20, 5291.10, 2329.27, 1112.14, 216.758, 94.3998, 8.71282]

# # ── Solve ─────────────────────────────────────────────────────────────────────
# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# # ── Plot ──────────────────────────────────────────────────────────────────────
# dem_safe  = np.clip(dem, 1e-40, None)
# log_dem   = np.log10(dem_safe)

# # Filter + intensity label for subtitle
# filter_label = ",  ".join(
#     f"{f} [{i:.4g}]" for f, i in zip(filters, intensities)
# )

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0,label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# # Titles
# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13,  pad=12,
# ) #fontweight="bold",
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# # Subtitle with filters + intensities (wrapped as a second text line)
# fig.text(
#     0.5, 0.91,
#     filter_label,
#     ha="center", va="top",
#     fontsize=8.5, color="black",
#     wrap=True,
# )

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# # ── Save ──────────────────────────────────────────────────────────────────────
# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")

# ##################################################################################################################################
# ##################################################################################################################################


# # import matplotlib.pyplot as plt
# # import numpy as np
# # from xrtpy.response.tools import generate_temperature_responses
# # from xrtpy.xrt_dem_iterative import XRTDEMIterative

# # ── DEM 15 metadata ──────────────────────────────────────────────────────────
# DEM_ID           = 15
# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"]
# intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174,    1.48818]

# # ── Solve ─────────────────────────────────────────────────────────────────────
# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# # ── Plot ──────────────────────────────────────────────────────────────────────
# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0,label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# # ── Save ──────────────────────────────────────────────────────────────────────
# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")


# ##### 15 with 6 filters

# # ── DEM 15 metadata ──────────────────────────────────────────────────────────
# DEM_ID           = 15
# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# filters     = ["Al-poly", "C-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [8228.90,   6206.10,  2633.15,   593.411,  243.889,  23.1174]

# # ── Solve ─────────────────────────────────────────────────────────────────────
# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq


# # ── Plot ──────────────────────────────────────────────────────────────────────
# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# # ── Save ──────────────────────────────────────────────────────────────────────
# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")




# ##################################################################################################################################
# ##################################################################################################################################
# ##################################      INDEX 65                              ###################################################
# ##################################################################################################################################
# ##################################################################################################################################
# ##################################################################################################################################


# # ── DEM 65 with 5 filters  metadata ──────────────────────────────────────────────────────────
# DEM_ID           = 65
# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [7020.42,   4331.36,   1007.78,  416.894,  40.2357]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT

# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])


# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")


# ##### INdex 65 with 6 filters

# DEM_ID           = 65
# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [10060.6,  7020.42,   4331.36,   1007.78,  416.894,  40.2357]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")



# ##################################################################################################################################
# ##################################################################################################################################
# ##################################      INDEX 99                              ###################################################
# ##################################################################################################################################
# ##################################################################################################################################
# ##################################################################################################################################


# # ── DEM 99 with  6filters  metadata ──────────────────────────────────────────────────────────

# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 99
# filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [87662.9,  59826.3,   35197.6,   8244.95,  3381.95,  338.260]


# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq


# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0,label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])


# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")

# ################################################################################
# ######################### Index 99 with 5 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 99
# filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [59826.3,   35197.6,   8244.95,  3381.95,  338.260]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")



# ################################################################################
# ######################### Index 99 with 7 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 99
# filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"]
# intensities = [87662.9,  59826.3,   35197.6,   8244.95,  3381.95,  338.260,    22.4599]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")






# ##################################################################################################################################
# ##################################################################################################################################
# ##################################      INDEX 112                              ###################################################
# ##################################################################################################################################
# ##################################################################################################################################
# ##################################################################################################################################


# # ── DEM 112 with  5 filters  metadata ──────────────────────────────────────────────────────────

# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 112
# filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick"]
# intensities = [215903.,  93715.4,  8685.49,   330.326,   8099.61]


# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq


# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])


# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")

# ################################################################################
# ######################### Index 112 with 6 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 112
# filters     = ["Al-thick", "Be-thick", "Al-poly/Al-mesh", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"]
# intensities = [8685.49,   330.326,   3295230.,           8099.61,            293.908,            1253340.]



# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")



# ################################################################################
# ######################### Index 112 with 7 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 112
# filters     = ["Be-med", "Al-med", "Al-thick", "Be-thick", "Al-poly/Al-thick", "Al-poly/Be-thick", "C-poly/Ti-poly"]
# intensities = [215903.,  93715.4,  8685.49,   330.326,   8099.61,            293.908,            1253340.]


# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")




# ##################################################################################################################################
# ##################################################################################################################################
# ##################################      INDEX 159                              ###################################################
# ##################################################################################################################################
# ##################################################################################################################################
# ##################################################################################################################################


# # ── DEM 159 with  5 filters  metadata ──────────────────────────────────────────────────────────

# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 159
# filters     = ["Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [8086.46,   5438.57,   1621.40,  654.170,  72.3524]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq


# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])


# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")

# ################################################################################
# ######################### Index 159 with 6 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 159
# filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick"]
# intensities = [11054.7,  8086.46,   5438.57,   1621.40,  654.170,  72.3524]



# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq


# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")



# ################################################################################
# ######################### Index 159 with 7 filters ##############################
# ################################################################################


# OBSERVATION_DATE = "2008-10-10T00:00:00"
# DATE_LABEL       = "10-Oct-2008 00:00:00"

# DEM_ID      = 159
# filters     = ["C-poly", "Ti-poly", "Be-thin", "Be-med", "Al-med", "Al-thick", "Be-thick"]
# intensities = [11054.7,  8086.46,   5438.57,   1621.40,  654.170,  72.3524,    6.14447]

# responses = generate_temperature_responses(filters, OBSERVATION_DATE)

# dem_solver = XRTDEMIterative(
#     observed_channel=filters,
#     observed_intensities=intensities,
#     temperature_responses=responses,
#     monte_carlo_runs=0,
# )
# dem_solver.solve()

# dem  = dem_solver.dem
# logT = dem_solver.logT
# chisq = dem_solver.chisq

# dem_safe     = np.clip(dem, 1e-40, None)
# log_dem      = np.log10(dem_safe)
# filter_label = ",  ".join(f"{f} [{i:.4g}]" for f, i in zip(filters, intensities))

# fig, ax = plt.subplots(figsize=(10, 6))

# ax.step(logT, log_dem, where="mid", color="black", linewidth=2.0, label=f"{DATE_LABEL}  |  x^2 = {chisq:.4f}")

# ax.set_title(
#     f"Hinode/XRT/XRTpy — Base DEM — DEM Index {DEM_ID}",
#     fontsize=13, pad=12,
# )
# ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
# ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)

# fig.text(0.5, 0.91, filter_label, ha="center", va="top", fontsize=8.5, color="black")

# ax.legend(fontsize=11)
# ax.grid(visible=True, alpha=0.3)
# ax.set_xlim(logT.min(), logT.max())

# finite = np.isfinite(log_dem)
# if np.any(finite):
#     y = log_dem[finite]
#     pad = 0.25 * (y.max() - y.min() + 1e-6)
#     ax.set_ylim(y.min() - pad, y.max() + pad)

# plt.tight_layout(rect=[0, 0, 1, 0.91])

# outfile = f"dem_idx{DEM_ID}_{DATE_LABEL.replace(' ', '_').replace(':', '').replace('-', '')}_{len(filters)}filters.png"
# fig.savefig(outfile, dpi=150)
# plt.show()
# print(f"Saved: {outfile}")


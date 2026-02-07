"""
Iterative DEM with Hinode/XRT
=============================

This example shows how to compute a differential emission measure (DEM)
from Hinode/XRT multi-filter observations using the XRTpy iterative solver.

Inputs:
- Filter names
- Observation date (for calibration/response computation)
- Observed intensities in DN/s/pix
"""

import numpy as np

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# %%
# 1) Choose filters + observation date
filters = ["Be-med","Be-thin","Al-poly","Al-poly/Ti-poly","Ti-poly","Al-thick"]
observation_date = "2007-07-10T10:58"

# %%
# 2) Provide observed intensities (DN/s/pix)
# Replace these with values measured from your data (per filter).
intensities = [
    1134.523,
    76.072,
    2433.566,
    152.297,
    304.208,
    3.837,
    ]


# %%
# 3) Generate temperature response functions for the selected filters
responses = generate_temperature_responses(filters, observation_date)

# %%
# 4) Run the DEM solver (base solution)
dem_solver = XRTDEMIterative(
    observed_channel=filters,
    observed_intensities=intensities,
    temperature_responses=responses,
)

dem_solver.solve()

# %%
# 5) Inspect results
dem = dem_solver.dem
logT = dem_solver.logT

print("logT grid:", logT)
print("DEM:", dem)

# %%
# 6) Plot DEM
dem_solver.plot_dem()

# %%
# 7) Monte Carlo uncertainty estimate (optional)
dem_solver_mc = XRTDEMIterative(
    observed_channel=filters,
    observed_intensities=intensities,
    temperature_responses=responses,
    monte_carlo_runs=100,
)

dem_solver_mc.solve()
dem_solver_mc.plot_dem_mc()
"""
Iterative DEM with Hinode/XRT
=============================

This example shows how to compute a differential emission measure (DEM)
from Hinode/XRT multi-filter observations using the XRTpy iterative solver.

Inputs
------
- Filter names
- Observation date (for temperature response computation)
- Observed intensities in DN/s/pix
- Intensity uncertainties in DN/s/pix (optional)
"""
from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# %%
# 1) Choose filters + observation date
filters = ["Al-mesh", "Al-poly","Be-thin","Be-thick","Al-med","Al-poly/Ti-poly"]
observation_date = "2021-07-20T16:04"

# %%
# 2) Provide observed intensities (DN/s/pix)
# Replace these with values measured from your data (per filter).
intensities = [
    362.382339, #Al-mesh
    181.407312, #Al-poly
    148.933070, #Be-thin
    0.002636,   #Be-thick
    1.459830,   #Al-med
    45.390125  #Al-poly/Ti-poly
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
    # You may optionally provide `intensity_uncertainties` in DN/s/pix.
    # If not provided, a default model is used: max(0.03 * intensity, 2 DN/s/pix)
    # intensity_uncertainties=[56.7, 3.8, 121.7, 7.6, 15.2, 0.2], # example values
)

dem_solver.solve()

# %%
# 5) Inspect results
dem_solver.summary()

# Access results directly for your own analysis:
dem = dem_solver.dem        # DEM(T) array [cm^-5 K^-1]
logT = dem_solver.logT      # log10 temperature grid [K]

# %%
# 6) Plot DEM
dem_solver.plot_dem()

# %%
# 7) Monte Carlo uncertainty estimate (optional)
# Re-run the solver with Monte Carlo to estimate DEM uncertainties.
# Each run perturbs the observed intensities by their uncertainties
# and re-solves, producing an ensemble of DEM curves.
dem_solver_mc = XRTDEMIterative(
    observed_channel=filters,
    observed_intensities=intensities,
    temperature_responses=responses,
    monte_carlo_runs=100,
)

dem_solver_mc.solve()
dem_solver_mc.plot_dem_mc()
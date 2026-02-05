.. _xrtpy-dem-overview:

===================================
Differential Emission Measure (DEM)
===================================

.. contents::
    :local:
    :depth: 2

Introduction
------------

This page describes the XRTpy iterative DEM solver and the inputs and outputs needed to run it.

The differential emission measure (DEM) describes how much plasma is present
in the solar corona as a function of temperature. It is a key diagnostic for
understanding coronal heating, solar flares, and the thermal structure of
active regions.

Hinode/XRT is well suited for DEM analysis because it observes the corona
through multiple broadband filters, each sensitive to different temperature
ranges. By combining these channels, we can infer a DEM(T) that reproduces 
the observed filter intensities when combined with the instrument temperature 
response functions.

The solver is iterative in that it repeatedly adjusts a parameterized
DEM to minimize the difference between observed and modeled intensities.

DEM in XRTpy
------------
XRTpy provides a Python implementation of the iterative spline fitting method
originally available in IDL as `xrt_dem_iterative2.pro <https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro>`__.
The core solver is implemented in :class:`xrtpy.xrt_dem_iterative.dem_solver.XRTDEMIterative`.

Conceptually, the solver:
    1. Builds a regular grid in log10(T) between user-specified bounds.
    2. Interpolates the filter temperature responses onto that grid.
    3. Represents log10(DEM) as a spline in log10(T).
    4. Uses least-squares fitting (via ``lmfit``) to adjust the spline values so that the modeled filter intensities best match the observed intensities.
    5. Optionally performs Monte Carlo runs by perturbing the observed intensities with their errors and re-solving the DEM many times to estimate uncertainties.

This approach mirrors the structure and behavior of the IDL routine while providing a modern, 
fully open-source implementation in Python that integrates naturally with the scientific Python ecosystem.

Required inputs
---------------
The DEM workflow requires three main input pieces:

1. Observed channels (filters)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Type: ``list`` of ``str``
* Description: Names of the filters used in the observation, for example ``"Al-mesh"`` or ``"Be-thin"``.

2. Temperature response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The DEM class requires temperature response functions for each filter, which
describe the instrument sensitivity as a function of temperature. These
responses can be generated for a chosen observation date using utilities
provided in ``xrtpy.response.tools``.

* Units: DN s\ :sup:`-1` pix\ :sup:`-1` cm\ :sup:`5`
* Description: Instrument response as a function of temperature for each filter, matching the order of the filters.
* Can be generated using :func:`xrtpy.response.tools.generate_temperature_responses`.



Example
^^^^^^^

.. code-block:: python

    from xrtpy.response.tools import generate_temperature_responses

    filters = ["Al-poly", "C-poly", "Ti-poly"]
    responses = generate_temperature_responses(
        filters,
        "2012-07-10T12:03:20",
        )


3. Observed intensities
~~~~~~~~~~~~~~~~~~~~~~~
* Type: array-like
* Units: DN/s (normalized per pixel)
* Description: Measured intensities in each filter channel.
* Length must match the number of filters.


Overview of the XRTDEMIterative API
-----------------------------------
The main entry point is :class:`xrtpy.xrt_dem_iterative.dem_solver.XRTDEMIterative`.



Solving a DEM
~~~~~~~~~~~~~
.. code-block:: python

    from xrtpy.xrt_dem_iterative import XRTDEMIterative

    dem_solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )

    # Solve for the DEM
    dem_solver.solve()  # returns the DEM array, also stored in dem_solver.dem

    # Plot the DEM
    dem_solver.plot_dem()

    # Access DEM results
    dem = dem_solver.dem
    logT = dem_solver.logT

    # DEM solution
    print(dem_solver.dem)

    # Temperature grid
    print(dem_solver.logT)



Enabling Monte Carlo error estimates
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
To estimate uncertainties, you can enable Monte Carlo iterations. The solver
will perturb the observed intensities by their errors and re-solve the DEM
for each realization.

.. code-block:: python

    N_mc = 100  # number of Monte Carlo runs

    dem_solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=N_mc,
    )

    dem_solver.solve()

    # Monte Carlo DEM Plot
    dem_solver.plot_dem_mc()   # base DEM plus Monte Carlo curves

The arrays ``dem_solver.mc_dem``, ``dem_solver.mc_chisq``,
``dem_solver.mc_base_obs``, and ``dem_solver.mc_mod_obs`` are then available
for custom analysis.

Comparison with IDL
-------------------
The Python solver is designed to closely follow the logic of the
SolarSoft/IDL routine `xrt_dem_iterative2.pro <https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro>`__:

* Uses a regular log10(T) grid.
* Represents log10(DEM) at a set of spline knots.
* Uses a least-squares algorithm to minimize the chi-square statistic.
* Supports Monte Carlo noise realizations for uncertainty estimation.

Small numerical differences between the Python and IDL implementations can arise due to:

* Different interpolation choices (for example, cubic splines from SciPy).
* Differences in optimization libraries (lmfit versus IDL MPFIT).
* Floating-point rounding and platform-specific details.

Within these limits, the Python implementation is intended to produce
results that are consistent with the IDL tool.


Mathematical background
-----------------------
This section provides a short description of the equations solved by the
XRTpy DEM solver. It is intended for orientation rather than as a full
mathematical derivation.

The DEM inversion problem is mathematically ill posed, meaning that multiple
thermal distributions can reproduce the same set of observations. For each
filter channel :math:`i`, the observed intensity :math:`I_i` is related to the
DEM by

.. math::

    I_i = \int DEM(T)\, R_i(T)\, dT

where :math:`R_i(T)` is the temperature response function of the filter and
:math:`DEM(T)` describes the amount of emitting plasma as a function of
temperature.

Because the number of temperature bins typically exceeds the number of observed
channels, the inversion does not have a unique solution. XRTpy therefore uses a
forward-fitting approach. The DEM is represented as a smooth function in
:math:`\log_{10}(T)` using a small number of spline knots. Model intensities are
computed on a discrete temperature grid as

.. math::

    I_i^{model} = \sum_j DEM(T_j)\, R_i(T_j)\, T_j\, \Delta(\ln T)

The spline values are adjusted to minimize the chi-square statistic

.. math::

    \chi^2 = \sum_i \left[
        \frac{I_i^{model} - I_i^{obs}}{\sigma_i}
    \right]^2

where :math:`\sigma_i` are the observational uncertainties. Smoothness of the
solution is enforced implicitly through the spline representation and the
limited number of knots.

When Monte Carlo error estimation is enabled, the observed intensities are
perturbed according to their uncertainties,

.. math::

    I_i^{(k)} = I_i^{obs} + \mathcal{N}(0, \sigma_i)

and the DEM is re-fit for each realization. The spread of the resulting DEM
curves provides an estimate of the uncertainty in :math:`DEM(T)`.


Extended example with options
-----------------------------
Below is an extended example showing additional constructor options.
The values shown match the current defaults and are written out for clarity.

.. code-block:: python

    from xrtpy.response.tools import generate_temperature_responses
    from xrtpy.xrt_dem_iterative import XRTDEMIterative

    filters = ["Al-poly", "Ti-poly", "Be-thin", "C-poly"]
    #Example intensities
    intensities = [520.0, 104.0, 901.0, 458.0]  # DN/s/pix
    observation_date = "2012-10-27T10:00:03"

    responses = generate_temperature_responses(
        filters,
        observation_date,
    )

    dem_solver = XRTDEMIterative(
        observed_channel=filters,              # Filter names
        observed_intensities=intensities,      # Observed intensity values
        temperature_responses=responses,       # Instrument responses

        # Optional configuration:
        intensity_errors=None,                 # Observed uncertainties - default: auto-estimated (3%)
        minimum_bound_temperature=5.5,         # Minimum log T (default: 5.5)
        maximum_bound_temperature=8.0,         # Maximum log T (default: 8.0)
        logarithmic_temperature_step_size=0.1, # Bin width in log T (default: 0.1)
        monte_carlo_runs=100,                  # Number of Monte Carlo runs (default: none)
        max_iterations=2000,                   # Solver max iterations (default: 2000)
        normalization_factor=1e21,             # Normalization scaling factor (default: 1e21)
    )

    dem_solver.solve()
    dem_solver.plot_dem_mc()

.. note::
    The values shown above correspond to the solver defaults and are written
    out here to illustrate which parameters can be tuned. You can adjust these
    to suit your specific analysis needs. This mirrors the flexibility of the
    IDL routine ``xrt_dem_iterative2.pro``.

.. Acknowledgement
.. ---------------
.. *Development of the DEM solver in XRTpy has been supported in part by 
.. a NASA Heliophysics Tools and Methods (HTM) program grant (ROSES-2025, 
.. element B.20). This effort reflects the ongoing transition of DEM 
.. capabilities from legacy IDL routines into modern, open-source Python 
.. tools for the solar physics community.*


References
----------
The following references provide background on the Hinode/XRT instrument and
its use in coronal diagnostics:

- Golub, L., et al. (2004), Solar Physics, 243, 63. :cite:p:`golub:2004`
- Weber, M. A., et al. (2004), Astrophysical Journal, 605, 528. :cite:p:`weber:2004`

Notes and warnings
------------------
The solver performs basic validation of user inputs and may emit warnings in
cases where observed intensities are non-physical (e.g., negative values) or
approach known instrument limits. These warnings do not stop execution but are
intended to help users assess the reliability of the results.

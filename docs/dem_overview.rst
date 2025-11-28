.. _xrtpy-dem-overview:

===================================
Differential Emission Measure (DEM)
===================================

.. contents::
    :local:
    :depth: 2

Introduction
------------

The differential emission measure (DEM) describes how much plasma is present
in the solar corona as a function of temperature. It is a key diagnostic for
understanding coronal heating, solar flares, and the thermal structure of
active regions.

Hinode/XRT is well suited for DEM analysis because it observes the corona
through multiple broadband filters, each sensitive to different temperature
ranges. By combining these channels, we can infer a temperature distribution
DEM(T) that explains the observed X-ray intensities.


DEM in XRTpy
------------
XRTpy provides a Python implementation of the iterative spline fitting method
originally available in IDL as `xrt_dem_iterative2.pro <http://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro>`_.
The core solver is implemented in :class:`xrtpy.xrt_dem_iterative.XRTDEMIterative`.


Conceptually, the solver:
    1. Builds a regular grid in log10(T) between user-specified bounds.
    2. Interpolates the filter temperature responses onto that grid.
    3. Represents log10(DEM) as a spline in log10(T).
    4. Uses least-squares fitting (via ``lmfit``) to adjust the spline values so that the modeled filter intensities match the observed intensities.
    5. Optionally performs Monte Carlo runs by perturbing the observed intensities with their errors and re-solving the DEM many times to estimate uncertainties.


Required inputs
---------------
The DEM workflow requires three main input pieces:

1. Observed channels (filters)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Type: ``str`` or ``list`` of ``str``
* Description: Names of the filters used in the observation, for example ``"Al-mesh"`` or ``"Be-thin"``.
* These must correspond to filters understood by XRTpy and must match the provided temperature responses one-to-one.

2. Observed intensities
~~~~~~~~~~~~~~~~~~~~~~~
* Type: array-like
* Units: DN/s (normalized per pixel)
* Description: Measured intensities in each filter channel.
* Length must match the number of filters.


3. Temperature response functions
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
* Type: ``list`` of :class:`xrtpy.response.TemperatureResponseFundamental`
* Units: DN s\ :sup:`-1` pix\ :sup:`-1` cm\ :sup:`5`
* Description: Instrument response as a function of temperature for each filter, matching the order of the filters.
* Can be generated using :func:`xrtpy.response.tools.generate_temperature_responses`.



Example
-------
A simple example with two filters:


.. code-block:: python
    
    from xrtpy.response.tools import generate_temperature_responses

    filters = ["Al-poly", "Ti-poly"]
    responses = generate_temperature_responses(
        filters,
        "2012-10-27T00:00:00",
        )


Overview of the XRTDEMIterative API
-----------------------------------
The main entry point is :class:`xrtpy.xrt_dem_iterative.XRTDEMIterative`.

Constructor
~~~~~~~~~~~
.. code-block:: python

    from xrtpy.xrt_dem_iterative import XRTDEMIterative

    dem_solver = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        intensity_errors=None,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=8.0,
        logarithmic_temperature_step_size=0.1,
        monte_carlo_runs=0,
        max_iterations=2000,
        normalization_factor=1e21,
    )

    # Solve for the DEM
    dem = dem_solver.solve()  # returns the DEM array, also stored in dem_solver.dem

    # Plot the DEM
    dem_solver.plot_dem()


Enabling Monte Carlo error estimates
------------------------------------
To estimate uncertainties, you can enable Monte Carlo iterations. The solver
will perturb the observed intensities by their errors and re-solve the DEM
for each realization.

.. code-block:: python

    N_mc = 50  # number of Monte Carlo runs

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
SolarSoft/IDL routine `xrt_dem_iterative2.pro <https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro>`_:

* Uses a regular log10(T) grid.
* Represents log10(DEM) at a set of spline knots.
* Uses a least-squares algorithm to minimize chi-square.
* Supports Monte Carlo noise realizations for uncertainty estimation.


Small numerical differences can arise due to:

* Different interpolation choices (for example, cubic splines from SciPy).
* Differences in optimization libraries (lmfit versus IDL MPFIT).
* Floating-point rounding and platform-specific details.

Within these limits, the Python implementation is intended to produce
results that are consistent with the IDL tool.



Mathematical background
-----------------------
The DEM inversion problem is ill posed. For each filter channel i,
the observed intensity :math:`I_i` is related to the DEM through:

.. math::

    I_i = \int DEM(T) \, R_i(T) \, dT

where :math:`R_i(T)` is the temperature response function for the filter,
and :math:`DEM(T)` is the unknown thermal distribution.

Since the number of temperature bins typically exceeds the number of
observed channels, the inversion does not have a unique solution. The
XRTpy solver uses a forward-fitting approach:

1. Assume a parametric form for log10(DEM(T)) using spline knots.
2. Compute model intensities:

    .. math::

        I_i^{model} = \sum_j DEM(T_j)\, R_i(T_j)\, T_j\, \Delta(\ln T)

3. Adjust the spline values to minimize:

    .. math::

        \chi^2 = \sum_i \left[
            \frac{I_i^{model} - I_i^{obs}}{\sigma_i}
        \right]^2

Here :math:`\sigma_i` are the observational uncertainties. Smoothness and
the low number of spline knots help regularize the solution.


Monte Carlo iterations perturb the observed intensities:

.. math::

    I_i^{(k)} = I_i^{obs} + \mathcal{N}(0, \sigma_i)

and re-fit the DEM for each realization k. The spread in the resulting DEM
curves provides an estimate of the uncertainty in DEM(T).


Extended example with options
-----------------------------
Below is an extended example showing more constructor options explicitly.
These values match current defaults but are written out here for clarity.

.. code-block:: python

    from xrtpy.response.tools import generate_temperature_responses
    from xrtpy.xrt_dem_iterative import XRTDEMIterative

    filters = ["Al-poly", "Ti-poly", "Be-thin", "C-poly"]
    intensities = [2500.0, 1800.0, 900.0, 450.0]  # DN/s
    observation_date="2012-10-27T00:00:00"

    responses = generate_temperature_responses(
        filters,
        observation_date,
    )

    dem_solver = XRTDEMIterative(
        observed_channel=filters,              # Filter names
        observed_intensities=intensities,      # Observed intensity values
        temperature_responses=responses,       # Instrument responses

        # Optional configuration:
        intensity_errors=None,                 # Obs. uncertainties - default: auto-estimated (3%)
        minimum_bound_temperature=5.5,         # Minimum log T (default: 5.5)
        maximum_bound_temperature=8.0,         # Maximum log T (default: 8.0)
        logarithmic_temperature_step_size=0.1, # Bin width in log T (default: 0.1)
        monte_carlo_runs=100,                  # # of Monte Carlo runs (default: none)
        max_iterations=2000,                   # Solver max iterations (default: 2000)
        normalization_factor=1e21,             # Normalization saling factor (default: 1e21)
    )

    dem_solver.solve()
    dem_solver.plot_dem_mc()

.. note::
    The values shown above correspond to existing defaults in the solver, 
    but they are written out here to illustrate what can be tuned.  
    You can adjust these to best suit your analysis needs.  
    This mirrors the flexibility of the IDL routine 
    ``xrt_dem_iterative2.pro``.

.. Acknowledgement
.. ---------------
.. *Development of the DEM solver in XRTpy has been supported in part by 
.. a NASA Heliophysics Tools and Methods (HTM) program grant (ROSES-2025, 
.. element B.20). This effort reflects the ongoing transition of DEM 
.. capabilities from legacy IDL routines into modern, open-source Python 
.. tools for the solar physics community.*

References
----------
- Golub, L., et al. (2004), *Solar Physics*, 243, 63. :cite:`golub:2004`
- Weber, M. A., et al. (2004), *ApJ*, 605, 528. :cite:p:`weber:2004`.

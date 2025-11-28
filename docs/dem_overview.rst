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
        abundance_model="hybrid"
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



Comparison with IDL
-------------------
The Python solver mirrors the SolarSoft/IDL routine 
`xrt_dem_iterative2.pro <https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro>`_.  

While results are consistent, minor differences can occur due to 
interpolation choices and optimization details.


Mathematical Note: Ill-posed Nature of DEM Inversion
----------------------------------------------------
The DEM problem is inherently an **ill-posed mathematical inversion**.  

Given observed intensities :math:`I_i` in channels *i*, and their 
temperature response functions :math:`R_i(T)`, the relationship is:

.. math::

    I_i = \int DEM(T) \, R_i(T) \, dT

Recovering :math:`DEM(T)` from a small set of broadband channels is 
not unique and is technically fraught with perils.  

XRTpy (like the original IDL routine ``xrt_dem_iterative2.pro``) employs a 
**forward-fitting approach**:
- A trial DEM is guessed.
- It is folded through :math:`R_i(T)` to produce "model" intensities.
- The DEM spline points are adjusted to minimize chi-square between model and observed values.

Because the number of temperature bins typically exceeds the number 
of observations, the solution is constrained by assumptions (e.g., 
spline smoothness).  

Uncertainties are estimated through **Monte Carlo iterations**, where 
observations are perturbed by their errors and re-fit. The resulting 
distribution of DEM solutions gives an estimate of confidence.



Example Extension
-----------------
In addition to the required inputs, you can provide optional parameters 
to fine-tune the DEM solution.  
The example below shows all options explicitly set.  

.. code-block:: python

    from xrtpy.xrt_dem_iterative import XRTDEMIterative

    dem_solver = XRTDEMIterative(
        observed_channel=filters,         # Filter names
        observed_intensities=intensities, # Observed values
        temperature_responses=responses,  # Instrument responses

        intensity_errors=errors,   # Obs. uncertainties (default: 3%)
        min_T=5.6,                 # Min log T (default: 5.5)
        max_T=7.8,                 # Max log T (default: 8.0)
        dT=0.05,                   # Bin width in log T (default: 0.1)
        min_error=1.5,             # Minimum error floor (default: 2 DN)
        relative_error=0.02,       # Fractional error scaling (default: 0.03)
        monte_carlo_runs=50,       # # of Monte Carlo runs (default: none)
        max_iterations=3000,       # Solver max iterations (default: 2000)
        solv_factor=1e17,          # Scaling factor (default: 1e21)
    )

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
- Golub, L., et al. (2004), *Solar Physics*, 243, 63. :cite:p:`golub:2004`
- Weber, M. A., et al. (2004), *ApJ*, 605, 528. :cite:p:`weber:2004`.


.. Next Steps
.. ----------
.. - Explore example notebooks in the `examples/` directory. Coming soon. 

.. .. code-block:: python

..     from xrtpy.response.tools import generate_temperature_responses
..     from xrtpy.xrt_dem_iterative import XRTDEMIterative

..     # Define filters and observed intensities
..     filters = ["Al-poly","C-poly/Ti-poly"]
..     intensities = [250.0, 180.0]  # DN/s

..     # Generate responses
..     responses = generate_temperature_responses(
..         filters, 
..         observation_date="2007-07-10", 
..         abundance_model="hybrid"
..         )

..     # Solve XRT DEM
..     dem_solver = XRTDEMIterative(
..         observed_channel=filters,
..         observed_intensities=intensities,
..         temperature_responses=responses,
..     )

..     dem_result = dem_solver.solve()

..     dem_result.plot()

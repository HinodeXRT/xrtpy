.. _xrtpy-glossary:

********
Glossary
********

.. glossary::
   :sorted:

   CCD
      Charge-Coupled Device — the XRT's imaging detector, where photons are converted into electrons and then digitized into data numbers (DN).

   Channel
      A specific configuration of the XRT instrument, defined by the combination of filters in both filter wheels. Each channel affects the telescope's sensitivity and spectral response.

   Chi-Square (:math:`\chi^2`)
      A goodness-of-fit statistic that measures how well the modeled filter intensities match the observed intensities. It is computed as:

      .. math::

         \chi^2 = \sum_i \left( \frac{I_i^{\mathrm{model}} - I_i^{\mathrm{obs}}}{\sigma_i} \right)^2

      where :math:`I_i^{\mathrm{model}}` is the intensity predicted by the DEM forward
      model for filter :math:`i`, :math:`I_i^{\mathrm{obs}}` is the observed intensity,
      and :math:`\sigma_i` is the intensity uncertainty. A lower :math:`\chi^2` indicates
      a better fit. In :class:`~xrtpy.xrt_dem_iterative.XRTDEMIterative`, the
      :math:`\chi^2` of the observations is also a useful proxy for how well the recovered
      DEM matches the true DEM.

   Deconvolution
      A numerical image processing technique used to correct for the blurring caused by the telescope's Point Spread Function (PSF), improving sharpness and visibility of fine structures.

   DEM
      Differential Emission Measure (DEM) — a function that describes how much plasma is present along the line of sight as a function of temperature. See :ref:`xrtpy-dem-overview` for a detailed overview of DEM theory,
      usage, and the solver provided in XRTpy.

   DEM Inversion
      The process of determining the temperature distribution of coronal plasma (the DEM) from a small number of filter intensities. Since more temperature bins are used than available filters, the problem is mathematically
      underconstrained ("ill posed"), so regularization and smoothing are required to obtain a stable, physical solution.


   DN
      Data Number (DN) — the digital value recorded by the CCD, representing the detected photon flux, usually in DN s\ :sup:`−1`\ .

   Effective Area
      A measure of the telescope's throughput for a given channel, combining geometric area with efficiencies of mirrors, filters, and detector — all affected by contamination and aging.

   Emission Measure
      A quantity proportional to the square of the electron density integrated along the line of sight, used to infer the amount of emitting plasma.

   FilterRatio Method
      A technique used to estimate plasma temperature and emission measure by comparing the intensity of X-ray images taken with two different XRT filters. Assumes an isothermal plasma and relies on accurate temperature response functions.

   Light Leak
      Refers to periodic increases in visible stray light contamination in XRT images, likely caused by a pinhole or defect in the entrance filter. This results in unwanted signal contributions, particularly affecting synoptic composites and thin filters. Users are encouraged to consult the `SolarSoft XRT Analysis Guide`_ (see Table 2.2) and use `xrtpy.image_correction.remove_lightleak` for correction.


   Contamination (related to the XRT)
      Refers to the gradual accumulation of material on the CCD and focal plane filters (FPFs), which reduces instrument throughput. This time-dependent degradation impacts effective area calculations and must be accounted for in data analysis. Refer to Section 2.5.3 *Contamination* in the `SolarSoft XRT Analysis Guide`_ for more information.


   Monte Carlo DEM
      A method for estimating DEM uncertainty through repeated calculations using randomly perturbed input data. Each DEM solution is computed by
      adding random noise, based on the intensity uncertainties, to the observed intensities and then repeating the DEM inversion. The spread
      of the resulting solutions estimates the uncertainty in the DEM at each temperature. In :class:`~xrtpy.xrt_dem_iterative.XRTDEMIterative`, the number of
      realizations is controlled by the ``monte_carlo_runs`` parameter.

   PSF
      Point Spread Function — describes the response of the telescope to a point source of light. In XRTpy, it is used in deconvolution routines to sharpen images.

   Response matrix
      A two-dimensional array containing the temperature response of each XRT
      filter, interpolated onto the solver's regular :math:`\log_{10}(T)`
      temperature grid. This matrix connects the DEM to the modeled filter
      intensities through the forward model:

      .. math::

         I_i^{\mathrm{model}}
         = \sum_j \mathrm{DEM}(T_j)\,R_i(T_j)\,T_j\,\Delta(\ln T)


   Solar Emission Spectra
      Emission spectra produced by solar plasma across a range of temperatures, calculated using spectral models such as CHIANTI. These spectra are used in temperature response and filter ratio methods.

   Spline Knots
      Selected temperature locations used to parameterize the DEM curve in
      :class:`~xrtpy.xrt_dem_iterative.XRTDEMIterative`. The DEM is represented
      by a cubic spline, with the DEM values at a small number of knot locations
      in :math:`\log_{10}(T)` space optimized during fitting. The number of knots
      is set to

      .. math::

         \min\left(N_{\mathrm{filters}} - 1,\, 7\right),

      which limits the number of fitted spline parameters relative to the
      available observational constraints and matches the behavior of the IDL
      routine ``xrt_dem_iterative2.pro``.

   Temperature Response
      The expected response of the instrument to isothermal plasma as a function of temperature, given in units of DN cm\ :sup:`5` s\ :sup:`−1` pix\ :sup:`−1` for each filter channel. Calculated using CHIANTI atomic models and user-defined abundances.

.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf

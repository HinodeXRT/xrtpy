.. _xrtpy-glossary:

********
Glossary
********

.. glossary::
   :sorted:
   
   CCD
      Charge-Coupled Device — the XRT’s imaging detector, where photons are converted into electrons and then digitized into data numbers (DN).

   Channel
      A specific configuration of the XRT instrument, defined by the combination of filters in both filter wheels. Each channel affects the telescope’s sensitivity and spectral response.

   Deconvolution
      A numerical image processing technique used to correct for the blurring caused by the telescope’s Point Spread Function (PSF), improving sharpness and visibility of fine structures.
   
   DN
      Data Number (DN) — the digital value recorded by the CCD, representing the detected photon flux, usually in DN s\ :sup:`−1`\ .
   
   Effective Area
      A measure of the telescope's throughput for a given channel, combining geometric area with efficiencies of mirrors, filters, and detector — all affected by contamination and aging.

   Emission Measure
      A quantity proportional to the square of the electron density integrated along the line of sight, used to infer the amount of emitting plasma.

   Light Leak
      Refers to periodic increases in visible stray light contamination in XRT images, likely caused by a pinhole or defect in the entrance filter. This results in unwanted signal contributions, particularly affecting synoptic composites and thin filters. Users are encouraged to consult the `SolarSoft XRT Analysis Guide`_ (see Table 2.2) and use `xrtpy.image_correction.remove_lightleak` for correction.


   Contamination (related to the XRT)
      Refers to the gradual accumulation of material on the CCD and focal plane filters (FPFs), which reduces instrument throughput. This time-dependent degradation impacts effective area calculations and must be accounted for in data analysis. Refer to Section 2.5.3 *Contamination* in the `SolarSoft XRT Analysis Guide`_ for more information.
   
   PSF
      Point Spread Function — describes the response of the telescope to a point source of light. In XRTpy, it is used in deconvolution routines to sharpen images.
   
   Solar Emission Spectra
      Emission spectra produced by solar plasma across a range of temperatures, calculated using spectral models such as CHIANTI. These are used in temperature response and filter ratio methods.

   Temperature Response
      The expected response of the instrument to isothermal plasma as a function of temperature, given in units of DN cm\ :sup:`5` s\ :sup:`−1` pix\ :sup:`−1` for each filter channel. Calculated using CHIANTI atomic models and user-defined abundances.

.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf

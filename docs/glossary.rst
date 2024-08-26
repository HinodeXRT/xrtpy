.. currentmodule:: xrtpy

.. _glossary:

********
Glossary
********

.. glossary::
   :sorted:

   Light Leak
      In the context of the X-Ray Telescope (XRT), a light leak refers to periodic increases in visible stray light contamination, which is often categorized as stray light.
      This phenomenon is primarily attributed to a potential pinhole in the telescope's entrance filter, allowing some visible light transmission even when the visible light shutter is closed.
      This defect has resulted in significant visible light contributions to X-ray images, especially noticeable in some of the thinner XRT X-ray filters.
      Users of XRT data, particularly those working with affected filters, are encouraged to consult the `SolarSoft XRT Analysis Guide`_ for detailed analysis and
      potential corrections of the stray light in their datasets (see Table 2.2 in the guide for more details).

   DN
      Data number (DN) per unit.

   Contamination (related to the XRT)
      Contamination refers to the accumulation of contaminating material on the XRT CCD and focal plane filters (FPFs), which results in a decrease in sensitivity.
      Refer to Section 2.5.3 `Contamination` in the `SolarSoft XRT Analysis Guide`_ for more information about the XRT contamination.

   Temperature response
      Temperature response refers to the instrument's temperature response function for a specific filter channel. Units measured in
      DN cm\ :sup:`5` s\ :sup:`−1` pix\ :sup:`-1`\ .

   Solar-Emission-Spectra
      Solar emission spectra are a collection of plasma emission spectra corresponding to different temperatures. They provide information about the emitted radiation at various temperatures.

.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf

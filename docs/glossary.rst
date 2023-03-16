.. currentmodule:: xrtpy

.. _glossary:

********
Glossary
********

.. glossary::
   :sorted:

   CCD
      Charge-Coupled Device (CCD) camera onboard the XRT instrument.
      The XRT uses a back-illuminated three-phase CCD with 13.5 µm pixel-size and 2048×2048 array. Refer to Section 3.4 `CCD Camera System` in the `SolarSoft XRT Analysis Guide`_ for more information.

   Temperature (units)
      Temperature :math:`T`, is a `~astropy.units.Quantity` in unit of 'K'(e.g., degrees kelvin).

   Temperature response
      Instrument temperature response function, for a filter-channel. Units are
      DN cm\ :sup:`5` s\ :sup:`−1` pix\ :sup:`-1`\ .

   Solar-emission-spectra
      A set of plasma emission spectra for a set of temperature.



.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf

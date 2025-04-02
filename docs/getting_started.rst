===============
Getting Started
===============

XRTpy is a Python package developed for the analysis of observations made by the X-Ray Telescope (XRT)
on board the Hinode spacecraft. This page is intended for new users of XRTpy. For background information about XRT,
please refer to the `SolarSoft XRT Analysis Guide`_.

XRTpy Objects
*************

XRTpy currently provides access to the following core classes:

- `Channel`
- `Effective Area`
- `Temperature Response`

It also includes functionality to:

- Derive temperature and emission measure from image pairs,
- Sharpen images using the Point Spread Function (PSF),
- Correct synoptic images for the light leak (visible stray light).

Visit our Examples page for step-by-step Jupyter notebook guides on how to use each feature.

Channel
-------

The `Channel` class describes the configuration of a specific XRT filter channel. It includes details for the Charge-Coupled Device (CCD),
Entrance Filter, Focal Plane Filters, Geometry, and Mirrors.

Effective Area
--------------

XRTpy calculates the effective area for each XRT filter channel, accounting for time-dependent contamination on the CCD. For more details,
refer to the `SolarSoft XRT Analysis Guide`_.

Temperature Response
--------------------

XRTpy calculates the temperature response of XRT filter channels using the CHIANTI_ atomic database (version 10.0) and coronal abundances
(:cite:t:`feldman:1992`). This produces a response function as a function of temperature, using an assumed emission model
(:cite:t:`narukage:2011`, :cite:t:`narukage:2014`).

Deriving Temperature and Emission Measure
-----------------------------------------

The `temperature_from_filter_ratio` function allows you to derive plasma temperature and emission measure from a pair of XRT images using the filter-ratio method. This mirrors the logic in the SolarSoft IDL routine of the same name. A usage example is available in the Examples section.

Image Deconvolution with the PSF
--------------------------------

The `deconvolve` function applies image deconvolution using the instrument’s Point Spread Function (PSF) to sharpen XRT images. This is especially useful for recovering detail around bright or sharp solar structures.

Light Leak Correction
---------------------

The `remove_lightleak` function subtracts visible stray light from XRT synoptic composite images. This correction improves the quality of long-term coronal evolution studies. See our Examples section for how to use this function.

Abundance Model Options
-----------------------

By default, XRTpy uses CHIANTI coronal abundances (:cite:t:`feldman:1992`). You may also choose:

- `"hybrid"` — based on :cite:t:`Fludra:1999`
- `"photospheric"` — based on :cite:t:`Grevesse:2007`

To use a different abundance model:

.. code-block:: python

   from xrtpy.response import TemperatureResponseFundamental

   TemperatureResponseFundamental(
       'Al-poly',
       '2022-07-04T23:43:12',
       abundance_model='hybrid'
   )

You may also pass the `abundance_model` keyword to `temperature_from_filter_ratio`.

.. note::
   In the future, XRTpy may support additional emission model libraries beyond CHIANTI.

Data Products
*************

XRT data products are available through the XRT website. These include:

- `Level 1 Data`_ — Calibrated data using `xrt_prep`, in instrumental DN units.
- `Level 2 Synoptics`_ — Composite images from the synoptic observing program.

For more information, visit the `XRT data products`_ page.

X-Ray Filter Channels
*********************

XRT uses two filter wheels to configure the imaging filter channel. Each wheel includes several filters and an open slot:

Filter Wheel 1:
  - *Open*
  - Al-poly
  - C-poly
  - Be-thin
  - Be-med
  - Al-med

Filter Wheel 2:
  - *Open*
  - Al-mesh
  - Ti-poly
  - G-band (visible)
  - Al-thick
  - Be-thick

.. note::
   - Each wheel has an *Open* slot used when the filter is in the opposite wheel.
   - XRTpy does not support G-band image processing or response calculations.

Filter names in XRTpy are passed as strings like `'Ti-poly'`.

References
==========

Velasquez, J., Murphy, N., Reeves, K. K., Slavin, J., Weber, M., & Barnes, W. (2024).
*XRTpy: A Hinode-X-Ray Telescope Python Package*. JOSS, 9(100), 6396.
https://doi.org/10.21105/joss.06396

.. _CHIANTI: https://www.chiantidatabase.org/chianti_database_history.html
.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _Level 1 Data: https://xrt.cfa.harvard.edu/level1/
.. _Level 2 Synoptics: https://xrt.cfa.harvard.edu/data_products/Level2_Synoptics/
.. _XRT data products: https://xrt.cfa.harvard.edu/data_products/index.php
.. _xrt_prep: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _XRT temperature response with other choice of abundances: http://solar.physics.montana.edu/takeda/xrt_response/xrt_resp.html

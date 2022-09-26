===============
Getting Started
===============

XRTpy is a Python package being developed for the analysis of observations made by the X-Ray Telescope (XRT)
on the board Hinode spacecraft. This page is intended for new users of `xrtpy`. For more background information about XRT please refer to the `SolarSoft XRT Analysis Guide`_.


XRTpy Objects:
**************
XRTpy currently offers *Channel*, *Effective Area*, and
*Temperature Response* classes. These classes allow users
to explore properties of the XRT. Visit our Example page for detail example guides on how to use the XRTpy classes.


Channel
-------
Channel is an instrument configuration class that explores properties of the XRT. Channel offers a detailed review of instruments for a chosen
filter channel including the Charge-Coupled Device (CCD), Entrance Filter, Focus-Filter(s), Geometry, and Mirror(s).


Effective Area
--------------
XRTpy produces the effective areas for a set of XRT filter channels paired with thicknesses of the CCD contamination layer.
Refer to the `SolarSoft XRT Analysis Guide`_ for more information about the instrumental spectral responses.


Temperature Response
--------------------
XRTpy produces the temperature response for each XRT filter channel, assuming a spectral emission model, refer to :cite:t:`narukage:2011` and :cite:t:`narukage:2014`.
XRTpy default emission model is from CHIANTI. This structure contains data and information about a plasma emission model, as a function of wavelength and temperature.
The default model assumes coronal abundances :cite:t:`feldman:1992`.

.. note::
   XRTpy has future plans to accept other plasma emission spectra models.


X-Ray Filter Channel
*********************
The XRT controls filter imaging using two sequentially positioned filter wheels. Refer to Section 3 in the `X-Ray Telescope Instrument Guide`
in the `SolarSoft XRT Analysis Guide`_ for more information about the XRT filters. The existing filters are structured as so:

#. Filter Configuration
    #. Filter position
        #. Filter Wheel 1:
            -  *Open*
            -  Aluminum Polyimide (*Al-poly*)
            -  Carbon Polyimide (*C-poly*)
            -  Beryllium Thin (*Be-thin*)
            -  Beryllium Medium (*Be-med*)
            -  Aluminum Medium (*Al-med*)
        #. Filter Wheel 2:
            -  *Open*
            -  Aluminum Mesh (*Al-mesh*)
            -  Titanium Polyimide (*Ti-poly*)
            -  *G-band*
            -  Aluminum Thick (*Al-thick*)
            -  Beryllium Thick (*Be-thick*)
    #. *Open*
        Each filter wheel has an empty position, named 'Open'. The open position is in place when a filter on the other filter wheel is being used.
    #. *G-band*
        The G-Band filter allows visible light into the telescope and onto the CCD. XRTpy does not
        calculate the effective area or the temperature response for the G-Band filter.

.. note::
    Filters are expressed by their abbreviation when used in XRTpy. For example, if we want to explore the filter channel
    that selects the titanium-on-polyimide filter, then the string would be ``'Ti-poly'``. The process is the same for all XRT
    filter channels.

.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _xrt-cfa-harvard: https://xrt.cfa.harvard.edu/index.php

===============
Getting Started
===============

XRTpy is a Python package being developed for the analysis of observations made by the X-Ray Telescope (XRT)
on the board Hinode spacecraft. This page is intended for new users of `xrtpy`. For more background information about XRT please refer to the `SolarSoft XRT Analysis Guide`_.


XRTpy Objects:
**************
XRTpy currently offers `~xrtpy.response.channel.Channel`, *Effective Area*, and
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
XRTpy produces the temperature response for each XRT filter channel, assuming a spectral emission model, refer to `Narukage et al. 2011`_ and `Narukage et al. 2014`_.
XRTpy default emission model from CHIANTI. This structure contains data and information about a plasma emission model, as a function of wavelength and temperature.
The default model assumes coronal abundances `Feldman (1992)`_.

.. note::
   XRTpy has future plans to accept other plasma emission spectra model.


X-Ray Filter Channel
*********************
The XRT controls filter imaging using two sequentially positioned filter wheels, reference Section 3 `X-Ray Telescope Instrument Guide`
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
        Visible light shutter position. Reference the Section 3.5 `XRT Mechanisms` in the `SolarSoft XRT Analysis Guide`_ for more
        information about 'Open' shutter position.
    #. *G-band*
        The G-Band filter allows visible light into the telescope and onto the CCD.

.. note::
    Filters are expressed by their abbreviation when used in XRTpy. For example, if we want to explore the filter channel
    that selects the titanium-on-polyimide filter, then the string would be 'Ti-poly'. The process is the same for all XRT
    filter channels.

.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _xrt-cfa-harvard: https://xrt.cfa.harvard.edu/index.php

.. _Feldman (1992): https://doi.org/10.1088/0031-8949/46/3/002

.. _Narukage et al. 2011: https://doi.org/10.1007/s11207-010-9685-2
.. _Narukage et al. 2014: https://doi.org/10.1007/s11207-013-0368-7

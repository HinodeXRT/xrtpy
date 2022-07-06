===============
Getting Started
===============

XRTpy is a Python package being developed for the analysis of observations made by the X-Ray Telescope (XRT)
on the Hinode spacecraft. This page is intended for new users of `XRTpy`.

XRTpy Objects:
**************
XRTpy currently offers *Channel*, *Effective Area*, and *Temperature Response* classes. These are classes allow users
to explore properties of the XRT. Visit our Example page for detail example guides on how to use XRTpy.


Channel
-------
Channel is an instrument configuration class that explores properties of the (XRT). Channel offers a detailed review of instruments for a chosen
filter channel e.g. Charge-Coupled Device (CCD), Entrance Filter, Filter(s), Geometry, and Mirror(s). An example guide can be found in our Example page.


Effective Area
--------------
XRTpy produces the effective areas for a set of XRT X-Ray filter channels paired with thicknesses of the CCD contamination layer.
Reference the `SolarSoft XRT Analysis Guide`_ for more information about the instrumental spectral responses.
An example of how to calculate the effective areas can be found in our Example page.


Temperature Response
--------------------
XRTpy produces the temperature response for each XRT X-Ray filter channel, assuming a spectral emission model, reference `Narukage et al. 2011`_ and `Narukage et al. 2014`_.
XRTpy default emission model from CHIANTI. This structure contains data and information about a plasma emission model, as a function of wavelength and temperature.
The default model assumes coronal abundances `Feldman (1992)`_. An example of how to calculate the temperature response can be found in our Example page.

.. note::
   XRTpy has future plans to accept other plasma emission spectra model.


XRTpy Variable Functionalities
******************************
    #. Filters
        #. Practicable Filters
            #. XRTpy objects workable filter: Aluminum medium, Aluminum mesh, Aluminum thick, Aluminum polyimide, Beryllium medium, Beryllium thick, Beryllium thin, Carbon polyimide, and Titanium polyimide,
        #. Naming
            #. Filters are expressed by their abbreviation when used in a XRTpy. For example, if we want to explore the channel filter that selects the titanium-on-polyimide filter, then the string would be 'Ti-poly'. The process is the same for all XRT practicable filters.
    #. Distinguish instances - Filter & Mirror
        #. Filter instance in Filter Wheel
            #. The XRT data is recorded through nine X-Ray filters using two filter wheels. We are able to explore detailed information of the chosen XRT channel filter using `channel.filter_#`, where '#' is expressing filter wheel 1 or 2. For example, if we are exploring Carbon-on-Polyimide in filter wheel 1, we will be exploring channel.filter_1.
        #. Mirror
            #. XRTpy offers the ability to inspect the first and second surface mirror. To distinguish the mirrors we use `channel_mirror_#`, where '#' is the first or second mirror surface.

X-Ray Filter Channels
*********************
The XRT controls filter imaging using two sequentially positioned filter wheels, Figure 3.1 in the `SolarSoft XRT Analysis Guide`_ shows the XRT filter wheels as viewed from the sun.
The existing filters are structures as so:

#. Filter Configuration
    #. Filter position
        #. Filter Wheel 1:
            -  Open
            -  Aluminum Polyimide
            -  Carbon Polyimide
            -  Beryllium Thin
            -  Beryllium Medium
            -  Aluminum Medium
        #. Filter Wheel 2:
            -  Open
            -  Aluminum Mesh
            -  Titanium Polyimide
            -  G-band
            -  Aluminum Thick
            -  Beryllium Thick
    #. Open
        Visible light shutter position. Reference the XRT mechanisms in the `SolarSoft XRT Analysis Guide`_ for more
        information about 'Open' shutter position.
    #. G-band
        The G-Band filter allows visible light into the telescope and onto the CCD.


.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _xrt-cfa-harvard: https://xrt.cfa.harvard.edu/index.php

.. _Feldman (1992): https://doi.org/10.1088/0031-8949/46/3/002

.. _Narukage et al. 2011: https://doi.org/10.1007/s11207-010-9685-2
.. _Narukage et al. 2014: https://doi.org/10.1007/s11207-013-0368-7

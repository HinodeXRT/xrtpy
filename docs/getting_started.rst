===============
Getting Started
===============

XRTpy is a Python package being developed for the analysis of observations made by the X-Ray Telescope (XRT)
on the Hinode spacecraft. This page is intended for new users of `XRTpy`.

XRTpy Objects:
**************

XRTpy currently offers *Channel*, *Effective Area*, and *Temperature Response*.
Exaplin each object (what they do), detials on how to use them (Filter, date,time), Referecnes, and reference notebook/example.

Channel
-------
Channel explores properties of the X-Ray Telescope (XRT). Channel offers a detailed review of the instrument's for a chosen
channel filter e.g. Charge-Coupled Device (CCD), Entrance Filter, Filter(s), Geometry, and Mirror(s). An example guide can be found in our Example page.

Effective Area
--------------
XRTpy produces the effective areas for a set of XRT x-ray channels paired with thicknesses of the CCD contamination layer.
Reference the `SolarSoft XRT Analysis Guide`_ for more information about the instrumental spectral responses.
An example of how to calculate the effective areas can be found in our Example page.

Temperature Response
--------------------
XRTpy produces the temperature response for each XRT X-ray channel, assuming a spectral emission model, reference `Narukage et al. 2011`_ and `Narukage et al. 2014`_.
XRTpy default emission model from CHIANTI. This structure contains data and information about a plasma emission model, as a function of wavelength and temperature.
The default model assumes coronal abundances `Feldman (1992)`_. An example of how to calculate the temperature response can be found in our Example page.

.. note::
   XRTpy has future plans to accept other plasma emission spectra model.

Filters
*******
The XRT controls filter imaging using two sequentially positioned filter wheels, Figure 3.1: in the  `SolarSoft XRT Analysis Guide`_ shows the XRT filter wheels as viewed from the sun.
The existing filters are structures as so:

#. Filters
    #. Filter position
        #. Filter Wheel 1 position:
            -  Open
            -  Aluminum Polyimide
            -  Carbon Polyimide
            -  Beryllium thin
            -  Beryllium medium
            -  Aluminum med
    #. Filter Wheel 2 position
        -  Open
        -  Aluminum mesh
        -  Titanium Polyimide
        -  Gband
        -  Aluminum thick
        -  Beryllium thick
    #. Open
        Visible light shutter position. Reference the XRT mechanisms section 3.5 in the `SolarSoft XRT Analysis Guide`_ for more
        information about 'Open' shutter position.
    #. Naming
        Filters are expressed by their abbreviation when used in a XRTpy object. For example, if we want to explore the
        channel filter that selects the titanium-on-polyimide filter, then the string would be 'Ti-poly'. The process is the same for all xrt filters.

#. Mirror
    #. Explain XRT mirrors - how to use mirror 1 and 2.

.. note::
   The X-Ray Telescope Software Guide not intended to be a guide to use XRTpy.


.. _SolarSoft XRT Analysis Guide: https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf
.. _xrt-cfa-harvard: https://xrt.cfa.harvard.edu/index.php

.. _Feldman (1992): https://doi.org/10.1088/0031-8949/46/3/002

.. _Narukage et al. 2011: https://doi.org/10.1007/s11207-010-9685-2
.. _Narukage et al. 2014: https://doi.org/10.1007/s11207-013-0368-7

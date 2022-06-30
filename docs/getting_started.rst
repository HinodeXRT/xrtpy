===============
Getting Started
===============

XRTpy is a Python package being developed for the analysis of observations made by the X-Ray Telescope (XRT)
on the Hinode spacecraft. This page is intended for new users of `xrtpy`.


XRTpy Objects:
**************

XRTpy currently offers *Channel*, *Effective Area*, and *Temperature Response*.
Exaplin each object (what they do), detials on how to use them (Filter, date,time), Referecnes, and reference notebook/example.

Channel
-------
Test Channel text.

Effective Area
--------------
Testing Effective Area text.

Temperature Response
--------------------
Temperature Response Referecnes `Narukage et al. 2011`_ and `Narukage et al. 2014`_.


Filters
*******
The XRT controls filter imaging using two sequentially positioned filter wheels, Figure 3.1: in the  `SolarSoft XRT Analysis Guide`_ shows the XRT filter wheels as viewed from the sun.
The existing filters are structures as so:

#. Filters
    #. Filter Wheel 1 position:
        - Open
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



.. _Narukage et al. 2011: https://doi.org/10.1007/s11207-010-9685-2
.. _Narukage et al. 2014: https://doi.org/10.1007/s11207-013-0368-7

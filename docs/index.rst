.. _xrtpy-documentation:

*******************
XRTpy Documentation
*******************

.. image:: _static/images/XRTpy_logo.png
   :alt: XRTpy logo
   :align: center
   :scale: 50%

Welcome to the documentation for **XRTpy** (version |version|) — a Python_ package for analyzing data from the `X-Ray Telescope`_ (XRT) :cite:p:`golub:2007` aboard the Hinode_ spacecraft :cite:p:`kosugi:2007`.

XRTpy is designed for solar physicists, astronomers, students, and anyone curious about the Sun's dynamic outer atmosphere.
Whether you're conducting research or just beginning to explore the world of X-ray solar imaging, XRTpy gives you the tools to study the hot plasma of the solar corona using real space-based observations.

.. note::

   **New in v0.6.0:** XRTpy now includes
   :class:`~xrtpy.xrt_dem_iterative.XRTDEMIterative`, a Python implementation
   of the IDL routine ``xrt_dem_iterative2.pro`` for computing Differential
   Emission Measures from Hinode/XRT observations. See :ref:`xrtpy-dem-overview`
   and the :ref:`changelog <xrtpy-changelog>` for details.


.. toctree::
   :maxdepth: 1
   :caption: Contents:

   about_xrt
   install
   getting_started
   dem_overview
   generated/gallery/index
   reference/index
   acknowledging_xrtpy
   bibliography
   glossary
   feedback_communication
   contributing
   code_of_conduct
   changelog/index

Published Work
--------------

The following paper describes the XRTpy package and its initial release- v0.4.0:

:cite:p:`velasquez:2024`


Get Involved
------------

XRTpy is an open-source project and welcomes contributions from the community.

* **Source code:** `github.com/HinodeXRT/xrtpy <https://github.com/HinodeXRT/xrtpy>`__
* **Report a bug or request a feature:** `GitHub Issues <https://github.com/HinodeXRT/xrtpy/issues>`__
* **Contribute:** See the :ref:`contributing guide <xrtpy-contributing>`
* **Questions or feedback:** See the :ref:`feedback page <xrtpy-feedback-communication>`

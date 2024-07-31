---
title: 'XRTpy: A Hinode-X-Ray Telescope Python Package'

tags:
  - Python
  - Astronomy
  - X-Ray Telescope
  - Solar Physics
  - Helio

authors:
  - name: Joy Velasquez
    orcid: 0009-0005-4804-0946
    equal-contrib: true
    affiliation: "1"
  - name: Nicholas A. Murphy
    orcid: 0000-0001-6628-8033
    affiliation: "1"
    corresponding: true
  - name: Katharine K. Reeves
    orcid: 0000-0002-6903-6832
    affiliation: "1"
    equal-contrib: true
  - name: Jonathan Slavin
    orcid: 0000-0002-7597-6935
    affiliation: "1"
    corresponding: true
  - name: Mark Weber
    orcid: 0000-0001-7098-7064
    affiliation: "1"
    corresponding: true
  - name: Will T. Barnes
    orcid: 0000-0001-9642-6089
    affiliation: "2,3"

affiliations:
 - name: Center for Astrophysics | Harvard-Smithsonian 60 Garden Street, Cambridge, MA, USA
   index : 1
 - name: American University, 4400 Massachusetts Avenue NW, Washington, DC 20016, USA
   index : 2
 - name : NASA Goddard Space Flight Center, 8800 Greenbelt Road, Greenbelt, MD 20771, USA
   index : 3

date: 15 December 2023
bibliography: paper.bib
---

# Summary

The XRTpy Python package is a specialized tool developed for the analysis of observations made by the X-Ray Telescope (XRT) [@Golub:2007] aboard the Hinode spacecraft [@Kosugi:2007]. Hinode is a joint mission involving space agencies from Japan, the United States, and Europe, launched with the primary aim of providing multi-wavelength data from the photosphere to the upper corona, enabling continuous observations of the Sun. Within this mission, the XRT instrument stands out as a remarkable piece of technology, capable of capturing high-resolution images of the solar corona's hottest material, diagnosing coronal temperatures from less than 1 MK to more than 10 MK.

# Statement of need

XRTpy is a Python package developed for the analysis of observations from the X-Ray Telescope (XRT) aboard the Hinode spacecraft.
It offers a comprehensive range of functionalities, including object-oriented representation of instrument configuration, effective area calculations, temperature response computation, light leak subtraction, image sharpening, estimation of electron temperature, and emission measure derivation.
These capabilities empower researchers to explore and analyze XRT data comprehensively, contributing to a deeper understanding of solar phenomena.

The official analysis routines for Hinode are scripted in the Interactive Data Language (IDL).
The [SolarSoft XRT Analysis Guide](https://xrt.cfa.harvard.edu/resources/documents/XAG/XAG.pdf) serves as the official software and instrument guide for XRT data analysis.
XRTpy has been carefully written to ensure the consistency and replication of results obtained from the official IDL routines as described in the SolarSoft XRT Analysis Guide.
Although currently XRTpy does not have all the capabilities of the SolarSoft routines for XRT, the package is in continual development and more functionality will be added in the future.
This alignment with established practices and standards aims to facilitate a seamless transition for researchers while harnessing the benefits of Python in solar data analysis.

A shift towards Python is underway within both NASA and the wider scientific community.
With XRTpy, Python users can efficiently analyze and process Hinode/XRT data, bridging the gap between traditional IDL routines and the increasing adoption of Python within the scientific community.
This transition not only enhances accessibility but also supports the broader trend in the scientific community toward Python-based data analysis tools, thereby fostering a collaborative and efficient environment for solar researchers.


# Package Structure

XRTpy is equipped with a range of capabilities tailored for the comprehensive analysis of XRT observations.
The package is structured into distinct modules, each serving a specific purpose:

 - `xrtpy.response.channel`: This module defines the `Channel` class which offers access to the properties of XRT filters, including information on the CCD, Entrance Filter, Focus-Filter(s), Geometry, and Mirror(s).

 - `xrtpy.response.effective_area`: XRTpy's capability to calculate effective areas for various XRT filter channels, combined with CCD contamination layer thickness information, is crucial for understanding instrumental spectral responses, as depicted in \autoref{fig:Figure 1}.

 - `xrtpy.response.temperature_response`: XRTpy provides the capability of computing the temperature response of all the XRT filter channels. It does this by relying on a spectral emission model, drawing from [@Narukage:2011] and [@Narukage:2014]. Users can choose from a range of CHIANTI abundance sets, including the default model with coronal abundances, hybrid, and photospheric options. \autoref{fig:Figure 2} shows XRTpy's temperature response calculations for all XRT filters across the different CHIANTI abundance sets. The CHIANTI database is described in [@Dere:1997] and the version used in XRTpy is CHIANTI version 10.0 [@Del_Zanna:2020]. Researchers have the flexibility to select the abundance model that best aligns with their research requirements.

 - `xrtpy.response.temperature_from_filter_ratio`: This module contains the `temperature_from_filter_ratio` function, which derives temperature and emission measure maps for a pair of images using the filter ratio method. \autoref{fig:Figure 3} illustrates an example usage of this function.

 - `xrtpy.image_correction.deconvolve`: Deconvolution is a powerful technique for improving image sharpness. The `deconvolve` function applies deconvolution to XRT images, effectively reducing blurring effects caused by the telescope's point spread function.

 - `xrtpy.image_correction.remove_lightleak`: The `remove_lightleak` function in this module eliminates light leak (visible stray light) from XRT synoptic composite images. This results in cleaner and more precise images suitable for in-depth analysis.

XRTpy supports multiple elemental abundance sets, including CHIANTI coronal abundances [@Feldman:1992], hybrid abundances (based on [@Fludra-and-Schmelz:1999] and [@Schmelz:2012]), and photospheric abundances (based on [@Grevesse:2007], [@Scott:2015], and [@Asplund:2009]).

XRTpy's capabilities are designed to empower researchers and scientists to fully exploit the potential of XRT data, offering the scientific community a unique opportunity to study the Sun's dynamic and complex behavior in a user-friendly and efficient manner.

![This graph displays the effective area for all X-ray focal-plane filters used in the XRT, plotted across a wavelength range of 0 to 70 angstroms. Each filter, represented by a unique color, shows distinct peaks that are important for choosing the best filter based on the wavelength being observed, and the curves demonstrate the instrument's ability to distinguish between different X-ray wavelengths. .\label{fig:Figure 1}](xrtpy_effective_area_plot.pdf)

![The temperature response (log scale) is plotted for all XRT X-ray focal-plane filters using XRTpy. Each curve represents the total instrument response as a function of temperature, integrated with different CHIANTI abundance models: Coronal (solid lines), Hybrid (dashed lines), and Photospheric (dotted lines). Highlights the sensitivity variations under different coronal conditions. .\label{fig:Figure 2}](xrtpy_temperature_response_plot.pdf)

![The application of the `temperature_from_filter_ratio` function is illustrated, demonstrating its role in calculating electron temperature and volume emission measure through filter ratios. The dataset, collected on January 28, 2011, between 01:31:55 and 01:32:05 UTC, comprises two images captured with specific filters. These images offer unique insights into solar conditions during the observed moments, as shown by [@Guidoni:2015]. .\label{fig:Figure 3}](xrtpy_temperature_from_filter_ratio_plot.pdf)

# Development of XRTpy Version 0.4.0

XRTpy version `0.4.0` was released on December 5, 2023. This version, available through the [Python Package Index](https://pypi.org/project/xrtpy/) (PyPI), can be installed using `pip` and is compatible with Python 3.9 and later.

In fostering collaboration within the solar physics community, interoperability with other packages is a work in progress that we've started on with the developers of [aiapy](https://aiapy.readthedocs.io/en/latest/) for SDO/AIA observations, [EISpack](https://eispac.readthedocs.io/en/latest/#) for Hinode-EUV imaging spectrometer (EIS) data analysis, and [irispy-lmsal](https://irispy-lmsal.readthedocs.io/en/stable/) for Interface Region Imaging Spectrograph (IRIS) observations. This integration provides users with a smooth and comprehensive analysis experience. Further building on the SunPy framework [@SunPy:2020], XRTpy effectively utilizes the `Map` object for handling Hinode/XRT image data.

The development of XRTpy is an open and collaborative effort, hosted on [Github](https://github.com/HinodeXRT/xrtpy) to ensure transparency and encourage community involvement. The project's documentation is comprehensive and continuously updated, available online via [Read the Docs](https://xrtpy.readthedocs.io/en/stable/). To maintain high-quality standards, XRTpy employs a robust testing framework built on [`pytest`](https://pytest.org) and [`GitHub Actions`](https://github.com/features/actions). This framework covers a range of aspects including different Python versions, online functionality, documentation integrity, software functionality, and code style checks, ensuring a reliable and effective tool for users.


# Acknowledgements

The development of XRTpy is supported by NASA contract NNM07AB07C to the Smithsonian Astrophysical Observatory.
Hinode is a Japanese mission developed and launched by ISAS/JAXA, with NAOJ as domestic partner and NASA and STFC (UK) as international partners.
It is operated by these agencies in co-operation with ESA and the NSC (Norway).
The XRTpy team would like to extend gratitude to the Python in Heliophysics Community at large for their contribution to XRTpy.
N.A.M. acknowledges support from NSF CSSI award 1931388.

# References

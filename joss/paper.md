---
title: 'XRTpy : A Hinode/X-Ray Telescope Python Package'

tags:
  - Python
  - Astronomy,
  - X-Ray Telescope
  - Solar Physics
  - Helio

authors:
  - name: Joy Velasquez
    orcid: 0009-0005-4804-0946
    equal-contrib: true
    affiliation: "1"
  - name: Nicholas Murphy
    orcid: 0000-0001-6628-8033
    affiliation: "1"
    corresponding: true
  - name: Katharine Reeves
    orcid: 0000-0002-6903-6832
    affiliation: "1"
    equal-contrib: true
  - name: Jonathan Slavin
    ordcid:0000-0002-7597-6935
    affiliation: "1"
    corresponding: true
  - name: Mark Weber
    orcid: 0000-0001-7098-7064
    affiliation: "1"
    corresponding: true

affiliations:
 - name: Center for Astrophysics | Harvard-Smithsonian 60 Garden Street. Cambridge, MA, USA
   index: 1

date: 12 December 2023
bibliography: paper.bib
---

# Summary

The XRTpy Python package is a specialized tool developed for the analysis of observations made by the X-Ray Telescope (XRT) [@Golub:2007] aboard the Hinode spacecraft [@Kosugi:2007]. Hinode, a joint mission involving space agencies from Japan, the United States, Europe, and the United Kingdom, was launched with the primary aim of providing multi-wavelength data from the photosphere to the upper corona, enabling continuous observations of the Sun. Within this mission, the XRT instrument stands out as a remarkable piece of technology, capable of capturing high-resolution images of the solar corona's hottest material, spanning temperatures from 1,000,000 to 10,000,000 Kelvin.

# Statement of need

XRTpy is a Python package developed for the analysis of observations from the X-Ray Telescope (XRT) aboard the Hinode spacecraft. It offers a comprehensive range of functionalities, including object-oriented representation of instrument configuration, effective area calculations, temperature response computation, light leak subtraction, image sharpening, electron temperature, and emission measure derivation, as well as various abundance model options. These capabilities empower researchers to explore and analyze XRT data comprehensively, contributing to a deeper understanding of solar phenomena.

The official analysis routines for Hinode are scripted in the Interactive Data Language (IDL), and the SolarSoft XRT Analysis Guide serves as the official software and instrument guide for XRT data analysis. XRTpy has been carefully written to ensure the consistency and replication of results obtained from the official IDL routines as described in the SolarSoft XRT Analysis Guide. Although currently XRTpy does not have all the capabilities of the XRT IDL routines, it is in continual development and will be adding more functionality in the future. This alignment with established practices and standards aims to facilitate a seamless transition for researchers while harnessing the benefits of Python in solar data analysis.

A shift towards Python is underway within both NASA and the wider scientific community. With XRTpy, Python users can efficiently analyze and process Hinode-XRT data, bridging the gap between traditional IDL routines and the increasing adoption of Python within the scientific community. This transition not only enhances accessibility but also supports the broader trend in the scientific community toward Python-based data analysis tools, thereby fostering a collaborative and efficient environment for solar researchers.


# Statement of Need

XRTpy is equipped with a range of capabilities tailored for the comprehensive analysis of XRT (X-Ray Telescope) data. The package is structured into distinct subpackages, each serving a specific purpose:

 - Channel: This subpackage, managed by the 'Channel' class, offers access to the properties of XRT filters, covering aspects like the CCD, Entrance Filter, Focus-Filter(s), Geometry, and Mirror(s).

 - Effective Area: XRTpy's capability to calculate effective areas for various XRT filter channels, combined with CCD contamination layer thickness information, is crucial for understanding instrumental spectral responses, as depicted in \autoref{fig:Figure 1}.

 - Temperature Response: XRTpy provides the capability of computing the temperature response of all the XRT filter channels. It does this by relying on a spectral emission model, drawing from [@Narukage:2011] and [@Narukage:2014]. Users can choose from a range of CHIANTI abundance sets, including the default model with coronal abundances detailed by [@Feldman:2014], Hybrid, and Photospheric options. \autoref{fig:Figure 2} shows XRTpy's temperature response calculations for all XRT filters across the different CHIANTI  abundance sets.

 - Deriving Temperature and Emission Measure: With the 'temperature_from_filter_ratio' routine, XRTpy simplifies the process of deriving temperature and emission measure maps for a pair of images using the filter ratio method, \autoref{fig:Figure 3} illustrates an example of the function application.

 - Enhancing Image Sharpness with Point Spread Function (Deconvolution): Deconvolution is a powerful technique for improving image sharpness. This tool applies deconvolution to XRT images, effectively reducing blurring effects caused by the telescope's point spread function.

 - Light Leak Correction: The 'remove_lightleak' function is designed to eliminate light leak (visible stray light) from XRT synoptic composite images. This results in cleaner and more precise images suitable for in-depth analysis.

 - Abundance Model: XRTpy supports multiple abundance sets, including CHIANTI Coronal abundances, Hybrid abundances (based on Fludra and Schmelz [1999]), and Photospheric abundances (based on [@Grevesse:2007]). Researchers have the flexibility to select the abundance sets model that best aligns with their research requirements.


![Figure 1: The Effective area for all XRT filters plotted using XRTpy. .\label{fig:Figure 1}](xrtpy_effective_area_plot.pdf)

![Figure 2: The temperature response is plotted for all XRT filters using XRTpy.  The plot also shows the effects of using different abundance models from Chianti for each filter.. .\label{fig:Figure 2}](xrtpy_temperature_response_plot.pdf)

![In Figure 3, the application of the 'temperature_from_filter_ratio' function is illustrated, demonstrating its role in calculating electron temperature and volume emission measure through filter ratios. The dataset, collected on January 28, 2011, between 01:31:55 and 01:32:05 UTC, comprises two images captured with specific filters. These images offer unique insights into solar conditions during the observed moments, as shown by [@Guidoni:2015]. .\label{fig:Figure 3}](xrtpy_temperature_from_filter_ratio_plot.pdf)

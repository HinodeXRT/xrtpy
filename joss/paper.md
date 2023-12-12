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

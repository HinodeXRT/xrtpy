.. _release guide:

**********************************************
Updating XRTpy testing data from SolarSoft IDL
**********************************************

This page is intentionally made for an XRTpy maintainer to update IDL testing text files to the XRTpy testing directory. The directory where the testing text files can be found is in xrtpy/xrtpy/response/test/data. 


IDL Testing Script 
==================



IDL testing data for channel
============================



IDL testing data for the effective area
=======================================
The IDL testing data for the effective area is located in test/data/effective_area_testing_files. This directory contains a folder for each corresponding filter channel e.g. Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly. Each filter folder contains text files. The text files contain header information i.e Filter , observation_data, and two rows of data; wavelength and effective area. The format of creating these text files is important for a proper testing run.  

title formatting 

Generating this data in SolarSoft IDL can be down running 

Mention how often we should test

Required change from the IDL create text file is to make sure all dashes are covered to underscore. 



IDL testing data for the temperature responses
==============================================
The IDL testing data for the temperature responses is located in test/data/temperature_responses_testing_files. This directory contains a folder for each corresponding filter channel e.g. Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly. Each filter folder contains text files. The text files contain header information i.e Filter , observation_data, and two rows of data; temperature and temperature responses. The format of creating these text files is important for a proper testing run.  

**********************************************
Updating XRTpy testing data from SolarSoft IDL
**********************************************

This page is intentionally made for an XRTpy maintainer to update testing text files from Solar Software Interactive Data Language (IDL) to the XRTpy testing environment for testing. A short explanation of the process consists of running IDL scripts with current dates to create new test files to test. Then moving those new testing files to the correct directories in xrtpy. Lastly,  running the testing code to test the new updated files. Consistent testing is crucial to xrtpy outputting the correct data for any requested date and time. 

IDL Testing Script 
==================
There is two IDL written testing scripts that will produce the testing text files. These IDL scripts are written in IDL format and meant to be ran in Solar Software Interactive Data Language (IDL). The IDL scripts are can be found at:  xrtpy/IDL_scritps/IDL_test_scripts. There is an effective area and temperature response script. These scripts contain information on the necessary changes to execute the code.

IDL testing data for channel
============================
After creating new testing text files from running the IDL scripts on your local machine. You will then transfer them to the correct directories in xrtpy. Following xrtpy/xrtpy/response/test/data you will find two additional folders effective_area_IDL_testing_files and temperature_response_IDL_testing_files. Place the testing files to their corresponding folder.

Moving up a directory xrtpy/xrtpy/response/test/data/effective_area_IDL_testing_files or temperature_response_IDL_testing_files you find additional folders for each filter.  Place the testing files to their corresponding filter folder. For each filter folder you will find text files with similar title naming. It is importing to maintain the same title naming conventions when updating new testing test files. Otherwise, the test code will not able to correctly find and test the update text files. Format should be (filter name first abbreviation)_(filter name second abbreviation)_YYYYMMDD_effective_area.txt. The text files contain header information i.e Filter , observation_data, and two rows of data; wavelength and effective area. The format of creating these text files is important for a proper test run. The IDL script will format this correctly.  

A short explanation of the process consists placing the new test text files in to the correct folder (effective area or temperature response ). Next, placing into the correct filter folder. Then making sure the title of the text file is in the correct format. Last, verifying the text files format is consistent. 

The purpose is to add new testing files to the already existing files.  This step can later be skipped by creating a code to automatically run monthly. 

IDL testing data for the effective area
=======================================
The IDL testing data for the effective area is located in response/test/data/effective_area_testing_files. This directory contains a folder for each corresponding filter channel e.g. Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly.  Each filter folder contains text files. The text files contain header information i.e Filter , observation_data, and two rows of data; temperature and  effective area.  Follow IDL Testing Scripts & IDL test files relocation to test the effective area. Once files are correctly placed. Run pytest test_effective_area.py to test. All results should pass for successful run. The purpose is to add new testing files to the already existing files.

IDL testing data for the temperature responses
==============================================
The IDL testing data for the temperature responses is located in test/data/temperature_responses_testing_files. This directory contains a folder for each corresponding filter channel e.g. Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly. Each filter folder contains text files. The text files contain header information i.e Filter , observation_data, and two rows of data; temperature and temperature responses. Follow IDL Testing Scripts & IDL test files relocation to test the effective area. Run pytest -v test_temperature_response.py to test. All results should pass for successful run. The purpose is to add new testing files to the already existing files.

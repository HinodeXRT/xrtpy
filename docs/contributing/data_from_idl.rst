Updating XRTpy Testing Data from SolarSoft IDL
==============================================

This page is designed for an XRTpy maintainer to update testing text files from the
Solar Software Interactive Data Language (IDL) to the XRTpy testing environment.
The process involves running IDL scripts with current dates to create new test files.
Subsequently, move these newly created testing files to the correct directories in XRTpy.
Finally, run the testing code to test the newly updated files. Consistent testing is
crucial for XRTpy to produce accurate data for any requested date and time.

IDL Testing Scripts
-------------------

There are two IDL testing scripts written to generate the testing text files. These IDL
scripts are written in IDL format and are intended to be run in Solar Software Interactive Data
Language (IDL). You can find the IDL scripts at `xrtpy/IDL_scripts/IDL_test_scripts`. There are
scripts for effective area and temperature response. These scripts contain information on the
necessary changes to execute the code.


IDL Testing Data for Channels
-----------------------------

After creating new testing text files by running the IDL scripts on your local machine, transfer them to the correct directories in XRTpy. Within `xrtpy/xrtpy/response/test/data`, you'll find two additional folders: `effective_area_IDL_testing_files` and `temperature_response_IDL_testing_files`. Place the testing files into their corresponding folders.

Moving up a directory to `xrtpy/xrtpy/response/test/data/effective_area_IDL_testing_files` or `temperature_response_IDL_testing_files`, you'll find additional folders for each filter. Place the testing files into their respective filter folders. Each filter folder contains text files with specific title formatting. It's important to maintain consistent title naming conventions when updating new testing files; the format should be (filter_name_first_abbreviation)_(filter_name_second_abbreviation)_YYYYMMDD_effective_area.txt. The text files contain header information such as Filter, observation_data, and two rows of data: wavelength and effective area. Properly formatting these text files is crucial for a successful test run. The IDL script will format this correctly.

The process involves placing the new test text files into the correct folder (effective area or temperature response), followed by placing them into the correct filter folder. Ensure the title of the text file is in the correct format and verify that the text file's format is consistent.

The purpose is to add new testing files to the existing ones. This step can later be automated by creating code to run monthly.


IDL Testing Data for Effective Area
-----------------------------------

The IDL testing data for the effective area is located in `response/test/data/effective_area_testing_files`. This directory contains a folder for each corresponding filter channel, e.g., Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly. Each filter folder contains text files with header information like Filter, observation_data, and two rows of data: temperature and effective area. Follow "IDL Testing Scripts & IDL Test Files Relocation" to test the effective area. Once files are correctly placed, run `pytest test_effective_area.py` to test. All results should pass for a successful run. The purpose is to add new testing files to the existing ones.

IDL Testing Data for Temperature Responses
------------------------------------------

The IDL testing data for temperature responses is located in `test/data/temperature_responses_testing_files`. This directory contains a folder for each corresponding filter channel, e.g., Al-mesh, Al-med, Al-poly, Al-thick, Be-med, Be-thick, Be-thin, C-poly, Ti-poly. Each filter folder contains text files with header information like Filter, observation_data, and two rows of data: temperature and temperature responses. Follow "IDL Testing Scripts & IDL Test Files Relocation" to test the temperature responses. Run `pytest -v test_temperature_response.py` to test. All results should pass for a successful run. The purpose is to add new testing files to the existing ones.

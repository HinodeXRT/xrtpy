import pytest

import numpy as np
import pkg_resources

import sunpy
import sunpy.map
import sunpy.io.special

from astropy import units as u
from datetime import datetime 
import glob

from xrtpy.response.temperature_response import TemperatureResponse


def get_IDL_data_files():
    
    directory = pkg_resources.resource_filename("xrtpy", 'response/tests/data/temperature_response_testing_files')
    
    filter_data_files = glob.glob(directory+'/**/*.txt') 

    return( sorted(filter_data_files) )

filenames = get_IDL_data_files()

def _IDL_raw_data_list(filename):

    with open(filename, "r") as filter_file:

        list_of_lists = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            list_of_lists.append(line_list)
    
    return list_of_lists

def IDL_test_filter_name(list_of_lists):
    return(str(list_of_lists[1][1]))

def IDL_test_date(list_of_lists):
    obs_date = str(list_of_lists[2][1])
    obs_time = str(list_of_lists[2][2])

    day = int(obs_date[0:2])
    
    month_datetime_object = datetime.strptime(obs_date[3:6], "%b")
    month = month_datetime_object.month

    year =int(obs_date[8:12])

    hour = int(obs_time[0:2])
    minute = int(obs_time[3:5])
    second = int(obs_time[6:8])

    observation_date_dt = datetime(year, month, day, hour, minute, second)

    return(observation_date_dt)


def _IDL_temperature_response_raw_data(filename):

    with open(filename, "r") as filter_file:

        list_of_lists = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            list_of_lists.append(line_list)
    
    new_list = []
    for list_ in list_of_lists[3:]:
        new_list.extend(list_)
    
    temperature_response = [float(i) for i in new_list]
    
    return(temperature_response)

@pytest.mark.parametrize("filename",filenames )

def test_temperature_response(filename,allclose):

    data_list = _IDL_raw_data_list(filename)
   
    filter_name = IDL_test_filter_name(data_list)
    filter_obs_date = IDL_test_date(data_list)

    IDL_temperature_response = _IDL_temperature_response_raw_data(filename)

    instance = TemperatureResponse(filter_name, filter_obs_date)
    actual_temperature_response = instance.temperature_response()

    assert allclose(actual_temperature_response.value , IDL_temperature_response,rtol=1e-6)

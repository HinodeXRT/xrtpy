import glob
import pkg_resources
import pytest

from datetime import datetime

from xrtpy.response.temperature_response import TemperatureResponseFundamental


def get_IDL_data_files():

    directory = pkg_resources.resource_filename(
        "xrtpy", "response/tests/data/temperature_response_IDL_testing_files"
    )

    filter_data_files = glob.glob(directory + "/**/*.txt")

    return sorted(filter_data_files)


filenames = get_IDL_data_files()


def _IDL_raw_data_list(filename):

    with open(filename, "r") as filter_file:

        IDL_data_list = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            IDL_data_list.append(line_list)

    return IDL_data_list


def IDL_test_filter_name(IDL_data_list):
    return str(IDL_data_list[1][1])


def IDL_test_date(IDL_data_list):
    obs_date = str(IDL_data_list[2][1])
    obs_time = str(IDL_data_list[2][2])

    day = int(obs_date[0:2])

    month_datetime_object = datetime.strptime(obs_date[3:6], "%b")
    month = month_datetime_object.month

    year = int(obs_date[8:12])

    hour = int(obs_time[0:2])
    minute = int(obs_time[3:5])
    second = int(obs_time[6:8])

    observation_date_dt = datetime(year, month, day, hour, minute, second)

    return observation_date_dt


def _IDL_temperature_response_raw_data(filename):

    with open(filename, "r") as filter_file:

        IDL_data_list = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            IDL_data_list.append(line_list)

    new_IDL_data_list = []
    for i in range(4, len(IDL_data_list)):
        new_IDL_data_list.append(IDL_data_list[i][1])

    temperature_response = [float(i) for i in new_IDL_data_list]

    return temperature_response


@pytest.mark.parametrize("filename", filenames)
def test_temperature_response(filename, allclose):

    IDL_data = _IDL_raw_data_list(filename)

    filter_name = IDL_test_filter_name(IDL_data)
    filter_obs_date = IDL_test_date(IDL_data)

    IDL_temperature_response = _IDL_temperature_response_raw_data(filename)

    instance = TemperatureResponseFundamental(filter_name, filter_obs_date)
    actual_temperature_response = instance.temperature_response()

    assert allclose(
        actual_temperature_response.value, IDL_temperature_response, rtol=1e-6
    )

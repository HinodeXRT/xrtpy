
import astropy.units as u
import numpy as np

import pytest


from astropy.utils.data import get_pkg_data_filenames
from datetime import datetime
from pathlib import Path

import pytest

from xrtpy.response.temperature_response import TemperatureResponseFundamental

_abundance_model_IDL_test_file_path = {
    "coronal": (
        Path(__file__).parent.absolute()
        / "data"
        / "temperature_response_coronal_IDL_testing_files"
    ),
    "hybrid": (
        Path(__file__).parent.absolute()
        / "data"
        / "temperature_response_hybrid_IDL_testing_files"
    ),
    "photospheric": (
        Path(__file__).parent.absolute()
        / "data"
        / "temperature_response_photospheric_IDL_testing_files"
    ),
}


def get_IDL_data_files():
    files = []
    data_root = "data/temperature_response_IDL_testing_files/"
    for top_dir in get_pkg_data_filenames(
        data_root,
        package="xrtpy.response.tests",
    ):
        files += list(
            get_pkg_data_filenames(
                top_dir, package="xrtpy.response.tests", pattern="*.txt"
            )
        )

    return sorted(files)


def get_IDL_data_files_older():
    path = (
        Path(__file__).parent.parent.absolute()
        / "tests"
        / "data"
        / "temperature_response_IDL_testing_files"
    )
    assert path.exists()
    filter_data_files = list(path.glob("**/*.txt"))
    # import pdb; pdb.set_trace()

    return sorted(filter_data_files)


filenames = get_IDL_data_files_older()

# filenames = get_IDL_data_files()


def _IDL_raw_data_list(filename):
    with open(filename) as filter_file:
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

    day = int(obs_date[:2])

    month_datetime_object = datetime.strptime(obs_date[3:6], "%b")
    month = month_datetime_object.month

    year = int(obs_date[8:12])

    hour = int(obs_time[:2])
    minute = int(obs_time[3:5])
    second = int(obs_time[6:8])
    return datetime(year, month, day, hour, minute, second)


def _IDL_temperature_response_raw_data(filename):
    with open(filename) as filter_file:
        IDL_data_list = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            IDL_data_list.append(line_list)

    new_IDL_data_list = [IDL_data_list[i][1] for i in range(4, len(IDL_data_list))]
    return [float(i) for i in new_IDL_data_list]


assert filenames


@pytest.mark.parametrize("filename", filenames)
def test_temperature_response(filename):  # allclose
    IDL_data = _IDL_raw_data_list(filename)

    filter_name = IDL_test_filter_name(IDL_data)
    filter_obs_date = IDL_test_date(IDL_data)

    IDL_temperature_response = _IDL_temperature_response_raw_data(filename)

    instance = TemperatureResponseFundamental(
        filter_name, filter_obs_date, abundance_model="coronal"
    )
    actual_temperature_response = instance.temperature_response()

    assert u.allclose(
        actual_temperature_response.value, IDL_temperature_response, rtol=1e-6
    ), filter_name

    # rtol=1e-1

    '''
    atol = actual_temperature_response.value.max() * 0.013
    assert allclose(
        actual_temperature_response.value,
        IDL_temperature_response,
        rtol=0.028,
        atol=atol,
    )
    '''


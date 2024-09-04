from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

from xrtpy.response.temperature_response import TemperatureResponseFundamental


def get_IDL_data_files(model):
    path = (
        Path(__file__).parent.parent.absolute()
        / "tests"
        / "data"
        / f"temperature_response_{model}_IDL_testing_files"
    )
    filter_data_files = list(path.glob("**/*.txt"))
    return sorted(filter_data_files)


def _IDL_raw_data_list(filename):
    with open(filename) as filter_file:  # noqa: PTH123
        IDL_data_list = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            IDL_data_list.append(line_list)
    return IDL_data_list


def IDL_test_abundance_name(IDL_data_list):
    return str(IDL_data_list[1][1])


def IDL_test_filter_name(IDL_data_list):
    return str(IDL_data_list[2][1])


def IDL_test_date(IDL_data_list):
    obs_date = str(IDL_data_list[3][1])
    obs_time = str(IDL_data_list[3][2])

    day = int(obs_date[:2])
    month_datetime_object = datetime.strptime(obs_date[3:6], "%b")
    month = month_datetime_object.month
    year = int(obs_date[8:12])

    hour = int(obs_time[:2])
    minute = int(obs_time[3:5])
    second = int(obs_time[6:8])
    return datetime(year, month, day, hour, minute, second)


def _IDL_temperature_response_raw_data(filename):
    with open(filename) as filter_file:  # noqa: PTH123
        IDL_data_list = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            IDL_data_list.append(line_list)

    new_IDL_data_list = [IDL_data_list[i][1] for i in range(5, len(IDL_data_list))]
    return [float(i) for i in new_IDL_data_list]


@pytest.mark.parametrize("abundance_model", ["coronal", "hybrid", "photospheric"])
def test_temperature_response(
    abundance_model,
):
    filenames = get_IDL_data_files(abundance_model)
    for filename in filenames:
        IDL_data = _IDL_raw_data_list(filename)
        filter_name = IDL_test_filter_name(IDL_data)
        filter_obs_date = IDL_test_date(IDL_data)
        IDL_temperature_response = _IDL_temperature_response_raw_data(filename)

        instance = TemperatureResponseFundamental(
            filter_name, filter_obs_date, abundance_model=abundance_model
        )

        actual_temperature_response = instance.temperature_response()
        atol = actual_temperature_response.value.max() * 0.013

        assert np.allclose(
            actual_temperature_response.value,
            IDL_temperature_response,
            rtol=0.028,
            atol=atol,
        )

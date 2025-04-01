from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filenames

from xrtpy.response.temperature_response import TemperatureResponseFundamental


def get_IDL_data_files(abundance):
    filter_data_files = []
    for dir in get_pkg_data_filenames(
        f"data/temperature_response_{abundance}_IDL_testing_files",
        package="xrtpy.response.tests",
    ):
        filter_data_files += list(Path(dir).glob("*.txt"))
    return sorted(filter_data_files)


filenames = (
    get_IDL_data_files("coronal")
    + get_IDL_data_files("hybrid")
    + get_IDL_data_files("photospheric")
)


@pytest.mark.parametrize("filename", filenames)
def test_temperature_response(filename):
    with Path.open(filename) as f:
        _ = f.readline()
        abundance = f.readline().split()[1]
        filter_name = f.readline().split()[1]
        filter_obs_date = " ".join(f.readline().split()[1:])

    # The IDL file might say "Sept" instead of "Sep"
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")

    IDL_data = np.loadtxt(filename, skiprows=5)
    IDL_temperature = IDL_data[:, 0] * u.K
    IDL_temperature_response = IDL_data[:, 1] * u.Unit("DN cm5 pix-1 s-1")

    instance = TemperatureResponseFundamental(
        filter_name,
        filter_obs_date,
        abundance_model=abundance,
    )

    actual_temperature_response = instance.temperature_response()

    IDL_temperature_response_interp = np.interp(
        instance.CHIANTI_temperature.value,
        IDL_temperature.value,
        IDL_temperature_response.value,
    ) * u.Unit("DN cm5 pix-1 s-1")

    assert u.allclose(
        actual_temperature_response,
        IDL_temperature_response_interp,
        rtol=0.005,
        atol=1e-6 * actual_temperature_response.max(),
    )

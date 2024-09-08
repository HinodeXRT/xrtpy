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
    # NOTE: Annoyingly the date strings use "Sept" instead of "Sep" for "September"
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    IDL_data = np.loadtxt(filename, skiprows=5)
    IDL_temperature = IDL_data[:, 0] * u.K
    IDL_temperature_response = IDL_data[:, 1] * u.Unit("DN cm5 pix-1 s-1")

    instance = TemperatureResponseFundamental(
        filter_name, filter_obs_date, abundance_model=abundance
    )

    IDL_temperature_response = np.interp(
        instance.CHIANTI_temperature, IDL_temperature, IDL_temperature_response
    )

    actual_temperature_response = instance.temperature_response()

    # NOTE: there may be small deviations where the response function is very small, likely
    # due to differences in the interpolation schemes. These are not critical as the response
    # is effectively zero in these regions anyway.
    i_valid = np.where(
        actual_temperature_response > 1e-8 * actual_temperature_response.max()
    )

    # NOTE: The relative tolerance is set comparatively high here because the CCD right gain
    # values in the IDL and xrtpy codes are explicitly different. See https://github.com/HinodeXRT/xrtpy/pull/76.
    # If the IDL results files are corrected to account for this updated gain, then this
    # relative tolerance can be set significantly lower. Setting the gains to be the same,
    # nearly all cases match within less than 1%.
    assert u.allclose(
        actual_temperature_response[i_valid],
        IDL_temperature_response[i_valid],
        rtol=5e-2,
    )

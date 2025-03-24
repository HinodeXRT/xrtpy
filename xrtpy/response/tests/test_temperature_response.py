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
# filenames = (
#     get_IDL_data_files("coronal")
# )


@pytest.mark.parametrize("filename", filenames)
def test_temperature_response(filename):

    with Path.open(filename) as f:
        _ = f.readline()
        abundance = f.readline().split()[1]       # e.g. 'coronal'
        filter_name = f.readline().split()[1]     # e.g. 'Ti_poly'
        filter_obs_date = " ".join(f.readline().split()[1:])  # e.g. '22 Sep 2015'

    # The IDL file might say "Sept" instead of "Sep"
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")


    IDL_data = np.loadtxt(filename, skiprows=5)
    IDL_temperature = IDL_data[:, 0] * u.K
    IDL_temperature_response = IDL_data[:, 1] * u.Unit("DN cm5 pix-1 s-1")

    instance = TemperatureResponseFundamental(
        filter_name,        # e.g. 'Ti_poly'
        filter_obs_date,    # e.g. '22 Sep 2015'
        abundance_model=abundance,  # e.g. 'coronal'
    )
    actual_temperature_response = instance.temperature_response()  # XRTpy

    IDL_temperature_response_interp = np.interp(
        instance.CHIANTI_temperature.value,  # XRTpy's temperature array
        IDL_temperature.value,               # IDL's temperature array
        IDL_temperature_response.value       # IDL's T-response array
    ) * u.Unit("DN cm5 pix-1 s-1")


    diff = np.abs((actual_temperature_response - IDL_temperature_response_interp)/ IDL_temperature_response_interp)
    
    max_diff = diff.max().value  # .value removes units for printing
    # Print if bigger than some threshold:
    if max_diff > 0.2:  # 20%
        print(f"[DEBUG] {filename} - max rel diff = {max_diff:.2%}")


    i_valid = np.where(actual_temperature_response > 1e-8 * actual_temperature_response.max())


    rtol = 0.3  #30% rtol
    atol = 1e-5 * actual_temperature_response.max()


    assert u.allclose(
        actual_temperature_response[i_valid],
        IDL_temperature_response_interp[i_valid],
        rtol=rtol,
        atol=atol,
    )
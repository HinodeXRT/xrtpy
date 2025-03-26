from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_filenames

from xrtpy.response.temperature_response import TemperatureResponseFundamental


import matplotlib.pyplot as plt

def plot_xrtpy_vs_idl(
    filename,
    chianti_temp, 
    xrtpy_response,
    IDL_temperature,
    IDL_temperature_response,
    IDL_temperature_response_interp
):
    # We'll do one figure, two lines:
    plt.figure()

    # X-axis: log10 of temperature in K
    x_xrtpy = np.log10(chianti_temp.value)

    # Plot XRTpy
    plt.plot(x_xrtpy, xrtpy_response.value, label="XRTpy response")

    # Plot IDL (interpolated) with a dashed line
    plt.plot(x_xrtpy, IDL_temperature_response_interp.value,
            "--", label="IDL (interp)")

    # Optionally, you could also show the raw IDL data 
    # (un-interpolated) on the same plot. 
    # That requires computing log10(IDL_temperature.value):
    x_idl_raw = np.log10(IDL_temperature.value)
    plt.plot(x_idl_raw, IDL_temperature_response.value, 
            "o", markersize=4, label="IDL raw (points)")

    plt.xlabel("log10( Temperature [K] )")
    plt.ylabel("Temperature Response [DN cm^5 pix^-1 s^-1]")
    plt.title(f"Debug plot:\n{filename}")
    plt.legend()
    plt.grid()
    plt.show()
    
# def maybe_debug_plot(
#     filename,
#     xrtpy_temp,
#     xrtpy_resp,
#     idl_temp,
#     idl_resp_raw,
#     idl_resp_interp,):
#     fig, ax = plt.subplots()

    # # X-axis is log10 of Temperature
    # ax.plot(
    #     np.log10(xrtpy_temp.value),
    #     xrtpy_resp.value,
    #     label="XRTpy response",
    #     linewidth=2,
    # )

    # ax.plot(
    #     np.log10(xrtpy_temp.value),
    #     idl_resp_interp.value,
    #     "--",
    #     label="IDL (interp)",
    #     linewidth=2,
    # )

    # ax.plot(
    #     np.log10(idl_temp.value),
    #     idl_resp_raw.value,
    #     "o",
    #     label="IDL raw (points)",
    #     markersize=4,
    # )

    # ax.set_xlabel("log10( Temperature [K] )")
    # ax.set_ylabel(r"Temperature Response [DN cm$^5$ pix$^{-1}$ s$^{-1}$]")
    # ax.set_title(f"Debug plot:\n{filename}")
    # ax.legend()

    # # Save figure to the same folder as the input file, but with .png
    # out_png = Path(filename).with_suffix(".png")
    # plt.savefig(out_png, dpi=150)
    # plt.close(fig)


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
        abundance = f.readline().split()[1]  # e.g. 'coronal'
        filter_name = f.readline().split()[1]  # e.g. 'Ti_poly'
        filter_obs_date = " ".join(f.readline().split()[1:])  # e.g. '22 Sep 2015'

    # The IDL file might say "Sept" instead of "Sep"
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")

    IDL_data = np.loadtxt(filename, skiprows=5)
    IDL_temperature = IDL_data[:, 0] * u.K
    IDL_temperature_response = IDL_data[:, 1] * u.Unit("DN cm5 pix-1 s-1")

    instance = TemperatureResponseFundamental(
        filter_name,  # e.g. 'Ti_poly'
        filter_obs_date,  # e.g. '22 Sep 2015'
        abundance_model=abundance,  # e.g. 'coronal'
    )
    actual_temperature_response = instance.temperature_response()  # XRTpy

    IDL_temperature_response_interp = np.interp(
        instance.CHIANTI_temperature.value,  # XRTpy's temperature array
        IDL_temperature.value,  # IDL's temperature array
        IDL_temperature_response.value,  # IDL's T-response array
    ) * u.Unit("DN cm5 pix-1 s-1")


    ###################################################################
    # diff = np.abs(
    #     (actual_temperature_response - IDL_temperature_response_interp)
    #     / IDL_temperature_response_interp
    # )
    # max_diff = diff.max().value  # .value removes units for printing
    # # Print if bigger than some threshold:
    # if max_diff > 0.2:  # 20%
    #     print(f"[DEBUG] {filename} - max rel diff = {max_diff:.2%}")
    ##################################################################
    i_valid = np.where(
        actual_temperature_response > 1e-8 * actual_temperature_response.max()
    )
    

    # rtol = .025#5e-2#0.1 
    # atol = 1e-7 * actual_temperature_response.max()
    rtol = 0.03#0.0275
    atol = 1e-7 * actual_temperature_response.max()

    assert u.allclose(
        actual_temperature_response[i_valid],
        IDL_temperature_response_interp[i_valid],
        rtol=rtol,
        atol=atol,
    )
    
    # Check the relative difference only where the response is "significant"
    # assert u.allclose(
    #     xrtpy_resp[i_valid],
    #     idl_resp_interp[i_valid],
    #     rtol=rtol,
    #     atol=atol,
    # ), f"XRTpy vs IDL mismatch for {filename}"

    # Optionally do a debug-plot every time (or only when mismatch above threshold):
    # maybe_debug_plot(
    #     filename,
    #     xrtpy_temp,
    #     xrtpy_resp,
    #     IDL_temperature,
    #     IDL_temperature_response,
    #     idl_resp_interp,
    # )
    
    ##############################################################################
    # success = u.allclose(
    #     actual_temperature_response[i_valid],
    #     IDL_temperature_response_interp[i_valid],
    #     rtol=rtol, atol=atol
    # )

    # if not success:
    #     # If we want to see exactly what's going on, let's measure the max difference:
    #     diff = np.abs((actual_temperature_response - IDL_temperature_response_interp)
    #                 / IDL_temperature_response_interp)
    #     max_diff = diff.max().value  # dimensionless

    #     # If bigger than 2% difference, do a debug plot:
    #     if max_diff > 0.02:
    #         plot_xrtpy_vs_idl(
    #             filename,
    #             instance.CHIANTI_temperature, 
    #             actual_temperature_response,
    #             IDL_temperature,
    #             IDL_temperature_response,
    #             IDL_temperature_response_interp
    #         )
    #         # Possibly keep the plot open or let the test proceed. 
    #         # You can also raise an assertion here if you want the test to fail:
    #         raise AssertionError(f"Test failed: max diff = {max_diff*100:.2f}%")

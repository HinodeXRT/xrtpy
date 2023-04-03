import numpy as np
import pkg_resources
import pytest

from pathlib import Path
from sunpy.map import Map

from xrtpy.response.xrt_deconvolve import xrt_deconvolve

test_file = "L1_XRT20120605_215839.9.fits"
idl_result_file = "L1_XRT20120605_215839.9.deconv.fits"


def get_observed_data():
    directory = pkg_resources.resource_filename(
        "xrtpy", "response/tests/data/xrt_deconvolve_testing_files"
    )
    data_file = Path(directory) / test_file

    return data_file


def get_IDL_results_data():
    directory = pkg_resources.resource_filename(
        "xrtpy", "response/tests/data/xrt_deconvolve_testing_files"
    )
    results_file = Path(directory) / idl_result_file

    return results_file


def test_unbinned():
    """
    Test case where the image is at full resolution, i.e. chip_sum = 1
    """
    test_data = get_observed_data()
    in_map = Map(test_data)
    IDL_result_image = get_IDL_results_data()
    IDL_result = Map(IDL_result_image)
    out_map = xrt_deconvolve(in_map)
    assert np.allclose(out_map.data, IDL_result.data, atol=1.0e-7)

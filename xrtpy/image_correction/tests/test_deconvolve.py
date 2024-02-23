from pathlib import Path

import numpy as np
import pkg_resources
from sunpy.map import Map

from xrtpy.image_correction.deconvolve import deconvolve

test_file = "L1_XRT20120605_215839.9.fits"
idl_result_file = "L1_XRT20120605_215839.9.deconv.fits"
test_file_binned = "L1_XRT20210730_175810.1.fits"
idl_result_file_binned = "L1_XRT20210730_175810.1_deconv.fits"


def get_observed_data():
    directory = pkg_resources.resource_filename("xrtpy", "image_correction/tests/data")
    data_file = Path(directory) / test_file

    return data_file


def get_IDL_results_data():
    directory = pkg_resources.resource_filename("xrtpy", "image_correction/tests/data")
    results_file = Path(directory) / idl_result_file

    return results_file


def get_observed_binned_data():
    directory = pkg_resources.resource_filename("xrtpy", "image_correction/tests/data")
    data_file = Path(directory) / test_file_binned

    return data_file


def get_IDL_results_binned_data():
    directory = pkg_resources.resource_filename("xrtpy", "image_correction/tests/data")
    results_file = Path(directory) / idl_result_file_binned

    return results_file


def test_unbinned():
    """
    Test case where the image is at full resolution, i.e. chip_sum = 1
    """
    test_data = get_observed_data()
    in_map = Map(test_data)
    IDL_result_image = get_IDL_results_data()
    IDL_result = Map(IDL_result_image)
    out_map = deconvolve(in_map)
    assert np.allclose(out_map.data, IDL_result.data, atol=1.0e-7)


def test_binned():
    """
    Test case where the image is binned, i.e. has chip_sum = 2
    This case needs higher tolerances (atol, rtol) because the binned PSF for
    the IDL and python codes don't match
    """
    test_data = get_observed_binned_data()
    in_map = Map(test_data)
    IDL_result_image = get_IDL_results_binned_data()
    IDL_result = Map(IDL_result_image)
    out_map = deconvolve(in_map)
    assert np.allclose(out_map.data, IDL_result.data, atol=4.1, rtol=0.008)

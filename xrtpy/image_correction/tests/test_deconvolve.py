from pathlib import Path

import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_path
from sunpy.map import Map

from xrtpy.image_correction.deconvolve import deconvolve


@pytest.fixture()
def data_dir():
    return Path(get_pkg_data_path("data", package="xrtpy.image_correction.tests"))


@pytest.mark.parametrize(
    ("observed_file", "idl_file", "atol", "rtol"),
    [
        (
            "L1_XRT20120605_215839.9.fits",
            "L1_XRT20120605_215839.9.deconv.fits",
            0,
            1e-7,
        ),
        # Case where image is binned, i.e. has chip_sum = 2
        # Needs higher tolerances because the binned PSF for the IDL and python codes don't match
        (
            "L1_XRT20210730_175810.1.fits",
            "L1_XRT20210730_175810.1_deconv.fits",
            4.1,
            0.008,
        ),
    ],
)
def test_unbinned(data_dir, observed_file, idl_file, rtol, atol):
    """
    Test deconvolution against IDL results for binned and unbinned cases
    """
    in_map = Map(data_dir / observed_file)
    IDL_result = Map(data_dir / idl_file)
    out_map = deconvolve(in_map)
    assert np.allclose(out_map.data, IDL_result.data, rtol=rtol, atol=atol)

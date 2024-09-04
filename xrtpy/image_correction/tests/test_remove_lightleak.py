from pathlib import Path

import numpy as np
import pytest
from astropy.utils.data import get_pkg_data_path
from sunpy.map import Map

from xrtpy.image_correction.remove_lightleak import remove_lightleak

data_dir = Path(
    get_pkg_data_path(
        "data/light_leak_testing_data_files", package="xrtpy.image_correction.tests"
    )
)
composite_filenames = sorted(
    (data_dir / "xrtpy_lightleak_data_test_files").glob("comp_XRT*.fits")
)
IDL_filenames = sorted(
    (data_dir / "IDL_lightleak_corrected_data_test_files").glob("ll_comp_XRT*.fits")
)


@pytest.mark.parametrize(
    ("idlfile", "compfile"), list(zip(IDL_filenames, composite_filenames))
)
def test_lightleak(idlfile, compfile):
    IDL_map = Map(idlfile)
    input_map = Map(compfile)
    ll_removed_map_xrtpy = remove_lightleak(input_map)
    # Because of rebinning for full resolution images, the match is worse
    # between IDL created images and XRTpy ones. IDL's method of rebinning
    # is different from that used by sunpy.
    atol = 0.75 if input_map.data.shape == (2048, 2048) else 1e-5
    assert np.allclose(ll_removed_map_xrtpy.data, IDL_map.data, atol=atol)

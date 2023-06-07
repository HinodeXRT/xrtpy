import pytest

from pathlib import Path
from sunpy.map import Map

from xrtpy.util.xrt_remove_lightleak import xrt_remove_lightleak


def get_IDL_data_file():
    """
    The XRT composite fits file that been lightleak corrected in
    IDL have been rename to begin with "ll" referring to lightleak.
    """
    directory = (
        Path(__file__).parent.absolute()
        / "data"
        / "light_leak_testing_data_files"
        / "IDL_lightleak_corrected_data_test_files"
    )
    data_files = directory.glob("ll_comp_XRT*.fits")
    return sorted(data_files)


IDL_filenames = get_IDL_data_file()


def get_composite_data_files():
    """
    The XRT composite fits file are no corrected in IDL.
    These files will be corrected using XRTpy.
    """
    directory = (
        Path(__file__).parent.absolute()
        / "data"
        / "light_leak_testing_data_files"
        / "xrtpy_lightleak_data_test_files"
    )
    data_file = directory.glob("comp_XRT*.fits")
    return sorted(data_file)


composite_filenames = get_composite_data_files()

# Using zip as an iterator to pair the data files together. Trouble-free method to use in pytest-parametrize
data_files = list(zip(IDL_filenames, composite_filenames))


@pytest.mark.parametrize(["idlfile", "compfile"], data_files)
def test_lightleak(idlfile, compfile, allclose):
    IDL_map = Map(idlfile)
    input_map = Map(compfile)

    ll_removed_map_xrtpy = xrt_remove_lightleak(input_map)

    assert allclose(ll_removed_map_xrtpy.data, IDL_map.data, atol=1e-5)

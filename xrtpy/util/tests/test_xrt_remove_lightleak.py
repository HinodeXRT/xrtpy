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


def get_observed_data():
    directory = Path(__file__).parent.absolute() / "data"
    data_file = list(directory.glob("comp_XRT*.fits"))
    return data_file[0]


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


def get_IDL_results_data(date_time):
    directory = Path(__file__).parent.absolute() / "data"
    # We give the IDL results data file have the same name as the input
    # composite image data file but with comp replaced with llfixed
    data_file = list(directory.glob(f"llfixed_XRT{date_time}.fits"))
    return data_file[0]


def test_one_case():
    input_data = get_observed_data()
    print(f"input_data = {input_data}")
    _, date_time = str(input_data).split("XRT")
    date_time, secdec, _ = date_time.split(".")
    date_time = ".".join([date_time, secdec])
    print(f"date_time = {date_time}")
    IDL_data = get_IDL_results_data(date_time)
    print(f"IDL_data = {IDL_data}")
    input_map = Map(input_data)
    IDL_map = Map(IDL_data)
    ll_removed_map = xrt_remove_lightleak(input_map)
    import numpy as np

    assert np.allclose(ll_removed_map.data, IDL_map.data)

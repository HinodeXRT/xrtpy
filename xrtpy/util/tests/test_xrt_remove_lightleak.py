import numpy as np

from pathlib import Path
from sunpy.map import Map

from xrtpy.util.xrt_remove_lightleak import xrt_remove_lightleak


def get_observed_data():
    directory = Path(__file__).parent.absolute() / "data"
    data_file = list(directory.glob("comp_XRT*.fits"))
    return data_file[0]


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
    assert np.allclose(ll_removed_map.data, IDL_map.data)

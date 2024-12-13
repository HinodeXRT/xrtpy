from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
#from astropy.utils.data import get_pkg_data_filenames

from xrtpy.response.channel import Channel
from xrtpy.response.effective_area import EffectiveAreaFundamental

channel_names = [
    "Al-mesh",
    "Al-poly",
    "C-poly",
    "Ti-poly",
    "Be-thin",
    "Be-med",
    "Al-med",
    "Al-thick",
    "Be-thick",
    "Al-poly/Al-mesh",
    "Al-poly/Ti-poly",
    "Al-poly/Al-thick",
    "Al-poly/Be-thick",
    "C-poly/Ti-poly",
]

channel_single_filter_names = [
    "Al-mesh",
    "Al-poly",
    "C-poly",
    "Ti-poly",
    "Be-thin",
    "Be-med",
    "Al-med",
    "Al-thick",
    "Be-thick",
]

valid_dates = [
    datetime(year=2006, month=9, day=25, hour=22, minute=1, second=1),
    datetime(year=2007, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2009, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2010, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2012, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2015, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2017, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2019, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2020, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2021, month=9, day=23, hour=22, minute=1, second=1),
    datetime(year=2022, month=9, day=23, hour=22, minute=1, second=1),
]

invalid_dates = [
    datetime(year=2006, month=8, day=25, hour=22, minute=1, second=1),
    datetime(year=2005, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2002, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=2000, month=9, day=22, hour=22, minute=1, second=1),
    datetime(year=1990, month=9, day=22, hour=22, minute=1, second=1),
]


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_name(channel_name):
    channel = Channel(channel_name)
    assert channel.name == channel_name


@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_filter_name(name):
    instance = EffectiveAreaFundamental(
        name, datetime(year=2013, month=9, day=22, hour=22, minute=0, second=0)
    )
    assert instance.name == name


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_contamination_on_CCD(name, date):
    instance = EffectiveAreaFundamental(name, date)
    assert 0 <= instance.contamination_on_CCD <= 1206


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_single_filter_names)
def test_EffectiveArea_contamination_on_filter(name, date):
    instance = EffectiveAreaFundamental(name, date)
    assert 0 <= instance.contamination_on_filter <= 2901


@pytest.mark.parametrize("date", invalid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_exception_is_raised(name, date):
    with pytest.raises(ValueError, match="Invalid date"):
        EffectiveAreaFundamental(name, date)


# def get_IDL_data_files(): 
#     filter_data_files = []
#     for dir in get_pkg_data_filenames(
#         "data/effective_area_IDL_testing_files", package="xrtpy.response.tests" 
#     ):
#         print(dir)
#         filter_data_files += list(Path(dir).glob("*.txt"))
#     return sorted(filter_data_files)

def get_IDL_data_files():
    data_dir = Path(__file__).parent / "data" / "effective_area_IDL_testing_files" 
    assert data_dir.exists(), f"Data directory {data_dir} does not exist."
    files = sorted(data_dir.glob("**/*.txt")) 
    #print(f"\n\n\nFound files: {files}\n\n")  # Debugging output
    return files



@pytest.mark.parametrize("filename", get_IDL_data_files())
def test_effective_area_compare_idl(filename):
    print(f"\n\nTesting file: {filename}\n")
    
    # Read the filter name and observation date from the file
    with filename.open() as f:
        filter_name = f.readline().split()[1]
        filter_obs_date = " ".join(f.readline().split()[1:])
        print(f"Filter name: {filter_name}, Observation date: {filter_obs_date}")  # Debugging output

    # Correct non-standard date format
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    
    # Load IDL data from the file
    IDL_data = np.loadtxt(filename, skiprows=3)
    IDL_wavelength = IDL_data[:, 0] * u.AA
    IDL_effective_area = IDL_data[:, 1] * u.cm**2
    #print(f"Loaded IDL data: {IDL_data.shape}\n")
    
    # Compute effective area using XRTpy
    instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
    actual_effective_area = instance.effective_area()
    #print(f"Calculated effective area (first 5): {actual_effective_area[:5]}\n")
    
    # Interpolate IDL effective area values to align with XRTpy wavelengths
    IDL_effective_area_interp = np.interp(instance.wavelength, IDL_wavelength, IDL_effective_area)
    #print(f"Interpolated IDL effective area (first 5): {IDL_effective_area_interp[:5]}\n")
    
    # Compare XRTpy and IDL effective areas using rtol=1e-4
    assert u.allclose(
        actual_effective_area,
        IDL_effective_area_interp,
        rtol=1e-4,
    ), f"Effective areas differ for filter {filter_name} on {filter_obs_date}"


## NOTE: This is marked as xfail because the IDL results that this test compares against
## are incorrect due to the use of quadratic interpolation in the contamination curves
## which leads to ringing near the edges in the contamination curve.
## See https://github.com/HinodeXRT/xrtpy/pull/284#issuecomment-2334503108
# @pytest.mark.xfail
# @pytest.mark.parametrize("filename", get_IDL_data_files())
# def test_effective_area_compare_idl(filename):
#     with Path.open(filename) as f:
#         filter_name = f.readline().split()[1]
#         filter_obs_date = " ".join(f.readline().split()[1:])
#     # NOTE: Annoyingly the date strings use "Sept" instead of "Sep" for "September"
#     filter_obs_date = filter_obs_date.replace("Sept", "Sep")
#     IDL_data = np.loadtxt(filename, skiprows=3)
#     IDL_wavelength = IDL_data[:, 0] * u.AA
#     IDL_effective_area = IDL_data[:, 1] * u.cm**2
#     instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
#     actual_effective_area = instance.effective_area()
#     IDL_effective_area = np.interp(
#         instance.wavelength, IDL_wavelength, IDL_effective_area
#     )
#     assert u.allclose(
#         actual_effective_area,
#         IDL_effective_area,
#         rtol=1e-4, #Moderate Precision
#     )
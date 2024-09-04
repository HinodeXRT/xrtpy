from datetime import datetime
from pathlib import Path

import numpy as np
import pytest
from astropy import units as u
from astropy.utils.data import get_pkg_data_filenames

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
    with pytest.raises(ValueError):  # noqa: PT011
        EffectiveAreaFundamental(name, date)


def get_IDL_data_files():
    filter_data_files = []
    for dir in get_pkg_data_filenames(
        "data/effective_area_IDL_testing_files", package="xrtpy.response.tests"
    ):
        filter_data_files += list(Path(dir).glob("*.txt"))
    return sorted(filter_data_files)


@pytest.mark.parametrize("filename", get_IDL_data_files())
def test_effective_area_compare_idl(filename):
    with Path.open(filename) as f:
        filter_name = f.readline().split()[1]
        filter_obs_date = " ".join(f.readline().split()[1:])
    # NOTE: Annoyingly the date strings use "Sept" instead of "Sep" for "September"
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")
    IDL_effective_area = np.loadtxt(filename, skiprows=3)[:, 1] * u.cm**2
    instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
    actual_effective_area = instance.effective_area()
    assert u.allclose(
        actual_effective_area,
        IDL_effective_area,
        atol=1e-2 * u.cm**2,
    )

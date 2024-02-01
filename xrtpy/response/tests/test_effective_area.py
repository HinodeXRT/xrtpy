from datetime import datetime
from pathlib import Path

import pytest
from astropy import units as u

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
    actual_attr_value = instance.name

    assert actual_attr_value == name


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
    with pytest.raises(ValueError):
        EffectiveAreaFundamental(name, date)


def get_IDL_data_files():
    directory = (
        Path(__file__).parent.parent.absolute()
        / "data"
        / "effective_area_IDL_testing_files"
    )
    filter_data_files = directory.glob("**/*.txt")
    return sorted(filter_data_files)


filenames = get_IDL_data_files()


def _IDL_raw_data_list(filename):
    with open(filename) as filter_file:
        list_of_IDL_effective_area_data = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            list_of_IDL_effective_area_data.append(line_list)

    return list_of_IDL_effective_area_data


def IDL_test_filter_name(list_of_lists):
    return str(list_of_lists[0][1])


def IDL_test_date(list_of_lists):
    obs_date = str(list_of_lists[1][1])
    obs_time = str(list_of_lists[1][2])

    day = int(obs_date[:2])

    month_datetime_object = datetime.strptime(obs_date[3:6], "%b")
    month = month_datetime_object.month

    year = int(obs_date[8:12])

    hour = int(obs_time[:2])
    minute = int(obs_time[3:5])
    second = int(obs_time[6:8])

    return datetime(year, month, day, hour, minute, second)


def _IDL_effective_area_raw_data(filename):
    with open(filename) as filter_file:
        list_of_lists = []
        for line in filter_file:
            stripped_line = line.strip()
            line_list = stripped_line.split()
            list_of_lists.append(line_list)

    effective_area = [list_of_lists[i][1] for i in range(3, len(list_of_lists))]
    effective_area = [float(i) for i in effective_area] * u.cm**2

    return effective_area


@pytest.mark.parametrize("filename", filenames)
def test_EffectiveAreaPreparatory_effective_area(filename, allclose):
    data_list = _IDL_raw_data_list(filename)

    filter_name = IDL_test_filter_name(data_list)
    filter_obs_date = IDL_test_date(data_list)

    IDL_effective_area = _IDL_effective_area_raw_data(filename)

    instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
    actual_effective_area = instance.effective_area()

    assert actual_effective_area.unit == IDL_effective_area.unit
    assert allclose(actual_effective_area.value, IDL_effective_area.value, atol=1e-2)

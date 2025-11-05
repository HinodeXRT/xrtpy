from datetime import datetime
from pathlib import Path

import numpy as np
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

# def test_ccd_contam_data_loads():
#     from xrtpy.response.effective_area import EffectiveAreaFundamental
#     eff = EffectiveAreaFundamental(...)  # use valid init args
#     data = eff.ccd_contam_data
#     assert isinstance(data, dict)  # or check for a specific key
# def test_ccd_contam_data_loads():
#     eff = EffectiveAreaFundamental(
#         "Al-poly", datetime(year=2013, month=9, day=22, hour=22, minute=1, second=1)
#     )
#     data = eff.ccd_contam_data
#     assert isinstance(data, dict)
#     assert "p1" in data and "p2" in data


def test_ccd_contam_data_loads():
    eff = EffectiveAreaFundamental(
        "Al-poly", datetime(year=2013, month=9, day=22, hour=22, minute=1, second=1)
    )
    data = eff.ccd_contam_data
    assert isinstance(data, dict)
    assert "p1" in data and "p2" in data


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
    # assert 0 <= instance.contamination_on_CCD <= 1206
    assert np.all(
        (instance.contamination_on_CCD >= 0) & (instance.contamination_on_CCD <= 1206)
    )


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_single_filter_names)
def test_EffectiveArea_contamination_on_filter1(name, date):
    instance = EffectiveAreaFundamental(name, date)
    # assert 0 <= instance.contamination_on_filter1 <= 2901
    assert np.all(
        (instance.contamination_on_filter1 >= 0)
        & (instance.contamination_on_filter1 <= 2901)
    )


@pytest.mark.parametrize("date", valid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_contamination_on_filter2(name, date):
    instance = EffectiveAreaFundamental(name, date)
    if instance.is_combo:
        result = instance.contamination_on_filter2
        assert result is None or np.all((result >= 0) & (result <= 2901))
        # assert result is None or (0 <= result <= 2901)
    else:
        assert instance.contamination_on_filter2 is None


@pytest.mark.parametrize("date", invalid_dates)
@pytest.mark.parametrize("name", channel_names)
def test_EffectiveArea_exception_is_raised(name, date):
    with pytest.raises(ValueError, match="Invalid date"):
        EffectiveAreaFundamental(name, date)


def get_IDL_data_files():
    data_dir = Path(__file__).parent / "data" / "effective_area_IDL_testing_files"
    assert data_dir.exists(), f"Data directory {data_dir} does not exist."
    files = sorted(data_dir.glob("**/*.txt"))
    # print(f"\n\n\nFound files: {files}\n\n")  # Debugging output
    return files


# NOTE: This is marked as xfail because the IDL results that this test compares against
# are incorrect due to the use of quadratic interpolation in the contamination curves
# which leads to ringing near the edges in the contamination curve.
# See https://github.com/HinodeXRT/xrtpy/pull/284#issuecomment-2334503108
# @pytest.mark.xfail
@pytest.mark.parametrize("filename", get_IDL_data_files())
def test_effective_area_compare_idl(filename):
    with filename.open() as f:
        filter_name = f.readline().split()[1]
        filter_obs_date = " ".join(f.readline().split()[1:])
    filter_obs_date = filter_obs_date.replace("Sept", "Sep")

    IDL_data = np.loadtxt(filename, skiprows=3)
    IDL_wavelength = IDL_data[:, 0] * u.AA
    IDL_effective_area = IDL_data[:, 1] * u.cm**2

    instance = EffectiveAreaFundamental(filter_name, filter_obs_date)
    actual_effective_area = instance.effective_area()

    XRTpy_effective_area = (
        np.interp(
            IDL_wavelength.value,
            instance.wavelength.value,
            actual_effective_area.value,
        )
        * u.cm**2
    )

    assert u.allclose(
        XRTpy_effective_area,
        IDL_effective_area,
        rtol=1e-4,
        atol=1.0e-2
        * u.cm**2,  # atol=1.0e-4 * u.cm**2, Revisit after getting Temp Resp working-JV
    )

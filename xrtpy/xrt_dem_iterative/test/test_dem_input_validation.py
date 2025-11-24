from importlib.resources import files
import numpy as np
import pytest
import sunpy
import sunpy.io.special
import sunpy.map


from xrtpy.response.channel import Channel
from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative


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


@pytest.mark.parametrize("channel_name", channel_names)
def test_channel_name(channel_name):
    channel = Channel(channel_name)
    assert channel.name == channel_name


# filename = Path(__file__).parent.parent.absolute() / "data" / "xrt_channels_v0017.genx"
filename = files("xrtpy.response.data") / "xrt_channels_v0017.genx"

v6_genx = sunpy.io.special.genx.read_genx(filename)
v6_genx_s = v6_genx["SAVEGEN0"]

_channel_name_to_index_mapping = {
    "Al-mesh": 0,
    "Al-poly": 1,
    "C-poly": 2,
    "Ti-poly": 3,
    "Be-thin": 4,
    "Be-med": 5,
    "Al-med": 6,
    "Al-thick": 7,
    "Be-thick": 8,
    "Al-poly/Al-mesh": 9,
    "Al-poly/Ti-poly": 10,
    "Al-poly/Al-thick": 11,
    "Al-poly/Be-thick": 12,
    "C-poly/Ti-poly": 13,
}


def test_dem_validate_inputs_basic():
    # Minimal “realistic” inputs for a DEM solve
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([2500.0, 1800.0])

    # Use real responses from XRTpy
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=8.0,
        logarithmic_temperature_step_size=0.1,
        monte_carlo_runs=0,
    )

    # Should NOT raise any error
    x.validate_inputs()


def test_dem_temperature_grid():
    filters = ["Al-poly"]
    intensities = np.array([1500.0])
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.1,
    )

    x.create_logT_grid()

    assert np.isclose(x.logT[0], 5.5)
    assert np.isclose(x.logT[-1], 7.5)
    assert len(x.logT) == 21  # (7.5-5.5)/0.1 + 1 = 21
    assert np.isclose(x.dlogT, 0.1)
    assert np.isclose(x.dlnT, np.log(10) * 0.1)

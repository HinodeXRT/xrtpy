from pathlib import Path

import numpy as np
import pytest
import sunpy
import sunpy.io.special
import sunpy.map

from xrtpy.response.channel import Channel

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


filename = Path(__file__).parent.parent.absolute() / "data" / "xrt_channels_v0017.genx"

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

def validate_inputs(self):
    """
    Run all internal validation checks again. Raises if any inputs are invalid.
    Useful for debugging or after programmatic changes.
    """
    # Check shape of intensity_errors
    if self._intensity_errors is not None:
        if self._intensity_errors.shape != self._observed_intensities.shape:
            raise ValueError("Length of intensity_errors must match observed_intensities.")

    # Check consistency between filters, intensities, and responses
    if not (
        len(self._observed_intensities)
        == len(self.responses)
        == len(self.observed_channel)
    ):
        raise ValueError(
            f"\nLength mismatch in inputs:\n"
            f"  Observed intensities: {len(self._observed_intensities)}\n"
            f"  Responses:            {len(self.responses)}\n"
            f"  Filter channels:      {len(self.observed_channel)}\n"
        )

    # Check temperature grid
    if self._dT <= 0:
        raise ValueError("dT must be a positive scalar.")

    for r in self.responses:
        logT_grid = np.log10(r.temperature.value)
        if not (self._min_T >= logT_grid.min() and self._max_T <= logT_grid.max()):
            raise ValueError(
                f"The specified temperature range [{self._min_T}, {self._max_T}] is outside the bounds of one or more filter response grids.\n"
                "Please ensure the temperature range fits within all responses."
            )


import pytest

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative


def test_validate_inputs_good_case():
    filters = ["Be-thin", "Be-med"]
    i_obs = [10000.0, 20000.0]
    resp = generate_temperature_responses(filters, obs_date="2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp)
    dem.validate_inputs()  # Should NOT raise

def test_validate_inputs_mismatched_errors():
    filters = ["Be-thin", "Be-med"]
    i_obs = [10000.0, 20000.0]
    i_err = [100.0]  # Wrong length
    resp = generate_temperature_responses(filters, obs_date="2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp, intensity_errors=i_err)
    with pytest.raises(ValueError, match="intensity_errors must match"):
        dem.validate_inputs()

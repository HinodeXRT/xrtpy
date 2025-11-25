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



# # -------------------------------------------------------------------------
# # 1) Valid configuration should pass validate_inputs
# # -------------------------------------------------------------------------

# def test_validate_inputs_valid_configuration(basic_responses, basic_intensities):
#     x = XRTDEMIterative(
#         observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#         observed_intensities=basic_intensities,
#         temperature_responses=basic_responses,
#         minimum_bound_temperature=5.5,
#         maximum_bound_temperature=8.0,
#         logarithmic_temperature_step_size=0.1,
#         monte_carlo_runs=0,
#         normalization_factor=1e21,
#     )

#     # Should not raise
#     x.validate_inputs()


# # -------------------------------------------------------------------------
# # 2) Empty observed_channel should raise
# # -------------------------------------------------------------------------

# def test_empty_observed_channel_raises(basic_responses, basic_intensities):
#     with pytest.raises(ValueError, match="`observed_channel` is required"):
#         XRTDEMIterative(
#             observed_channel=[],
#             observed_intensities=basic_intensities,
#             temperature_responses=basic_responses,
#         )


# # -------------------------------------------------------------------------
# # 3) Mismatched lengths of intensities / responses / channels
# # -------------------------------------------------------------------------

# def test_length_mismatch_raises():
#     responses = [DummyResponse("F1"), DummyResponse("F2")]
#     intensities = np.array([1000.0])  # only one value

#     with pytest.raises(ValueError, match="Length mismatch"):
#         XRTDEMIterative(
#             observed_channel=["F1", "F2"],
#             observed_intensities=intensities,
#             temperature_responses=responses,
#         )


# # -------------------------------------------------------------------------
# # 4) Temperature range outside response grid should raise
# # -------------------------------------------------------------------------

# def test_temperature_range_outside_responses_raises(basic_responses, basic_intensities):
#     # min T too low
#     with pytest.raises(ValueError, match="outside the bounds"):
#         XRTDEMIterative(
#             observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#             observed_intensities=basic_intensities,
#             temperature_responses=basic_responses,
#             minimum_bound_temperature=4.0,  # below dummy response range
#             maximum_bound_temperature=8.0,
#         )

#     # max T too high
#     with pytest.raises(ValueError, match="outside the bounds"):
#         XRTDEMIterative(
#             observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#             observed_intensities=basic_intensities,
#             temperature_responses=basic_responses,
#             minimum_bound_temperature=5.5,
#             maximum_bound_temperature=9.0,  # above dummy response range
#         )


# # -------------------------------------------------------------------------
# # 5) Negative or zero logarithmic_temperature_step_size should raise
# # -------------------------------------------------------------------------

# def test_negative_logarithmic_temperature_step_size_raises(basic_responses, basic_intensities):
#     with pytest.raises(ValueError, match="logarithmic_temperature_step_size must be a positive"):
#         XRTDEMIterative(
#             observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#             observed_intensities=basic_intensities,
#             temperature_responses=basic_responses,
#             minimum_bound_temperature=5.5,
#             maximum_bound_temperature=8.0,
#             logarithmic_temperature_step_size=-0.1,
#         )


# def test_too_few_temperature_bins_raises(basic_responses, basic_intensities):
#     # Choose a huge step so that fewer than 4 bins are produced
#     with pytest.raises(ValueError, match="Temperature grid must have at least 4 points"):
#         XRTDEMIterative(
#             observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#             observed_intensities=basic_intensities,
#             temperature_responses=basic_responses,
#             minimum_bound_temperature=5.5,
#             maximum_bound_temperature=5.8,
#             logarithmic_temperature_step_size=0.5,
#         )


# # -------------------------------------------------------------------------
# # 6) Monte Carlo runs validation
# # -------------------------------------------------------------------------

# def test_monte_carlo_runs_negative_raises():
#     with pytest.raises(ValueError, match="must be ≥ 0"):
#         make_iterative(monte_carlo_runs=-1)


# def test_monte_carlo_runs_bool_raises():
#     with pytest.raises(ValueError, match="must be a non-negative whole number, not a boolean"):
#         make_iterative(monte_carlo_runs=True)


# def test_monte_carlo_runs_float_non_integer_raises():
#     with pytest.raises(ValueError, match="Decimal values are not allowed"):
#         make_iterative(monte_carlo_runs=3.5)


# def test_monte_carlo_runs_zero_ok():
#     x = make_iterative(monte_carlo_runs=0)
#     assert x.monte_carlo_runs == 0


# def test_monte_carlo_runs_positive_integer_ok():
#     x = make_iterative(monte_carlo_runs=10)
#     assert x.monte_carlo_runs == 10


# # -------------------------------------------------------------------------
# # 7) Intensity errors validation in validate_inputs
# # -------------------------------------------------------------------------

# def test_intensity_errors_length_mismatch_raises(basic_responses, basic_intensities):
#     x = XRTDEMIterative(
#         observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#         observed_intensities=basic_intensities,
#         temperature_responses=basic_responses,
#         intensity_errors=np.array([1.0, 2.0]),  # wrong length
#     )

#     with pytest.raises(ValueError, match="Length of intensity_errors must match"):
#         x.validate_inputs()


# def test_intensity_errors_negative_raises(basic_responses, basic_intensities):
#     x = XRTDEMIterative(
#         observed_channel=["Filter-1", "Filter-2", "Filter-3"],
#         observed_intensities=basic_intensities,
#         temperature_responses=basic_responses,
#         intensity_errors=np.array([1.0, -2.0, 3.0]),
#     )

#     with pytest.raises(ValueError, match="must be finite and >= 0"):
#         x.validate_inputs()


# #Test to add later
# #both should be True
# # np.allclose(x.intensities_scaled,
# #             x.observed_intensities.value / x.normalization_factor)

# # np.allclose(x.sigma_scaled_intensity_errors,
# #             x.intensity_errors.to_value(u.DN/u.s) / x.normalization_factor)


from importlib.resources import files

import astropy.units as u
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
    resp = generate_temperature_responses(filters, "2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp)
    dem.validate_inputs()  # Should NOT raise


def test_validate_inputs_mismatched_errors():
    filters = ["Be-thin", "Be-med"]
    i_obs = [10000.0, 20000.0]
    i_err = [100.0]  # Wrong length
    resp = generate_temperature_responses(filters, "2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp, intensity_errors=i_err)
    with pytest.raises(ValueError, match="intensity_errors must match"):
        dem.validate_inputs()


def test_create_logT_grid():
    # Simple single-channel setup
    filters = ["Al-poly"]
    intensities = np.array([1500.0])
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    # Construct DEM object
    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.1,
    )

    # Create grid
    x.create_logT_grid()

    # 1 — Correct start and end
    assert x.logT[0] == pytest.approx(5.5)
    assert x.logT[-1] == pytest.approx(7.5)

    # 2 — Correct number of bins: (7.5 - 5.5)/0.1 + 1 = 21
    assert len(x.logT) == 21
    assert x.n_bins == 21

    # 3 — Correct spacing (uniform)
    diffs = np.diff(x.logT)
    assert np.allclose(diffs, 0.1, atol=1e-12)

    # 4 — dlogT and dlnT correct
    assert x.dlogT == pytest.approx(0.1)
    assert x.dlnT == pytest.approx(np.log(10) * 0.1)

    # 5 — T = 10**logT
    assert np.allclose(x.T.to_value(u.K), 10**x.logT)

    # 6 — logT strictly increasing
    assert np.all(np.diff(x.logT) > 0)


def test_estimate_initial_dem():

    # Step 1: Simple DEM case
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([1500.0, 2300.0])
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.1,
    )

    # Step 2: Create temperature grid & response matrix
    x.create_logT_grid()
    x._interpolate_responses_to_grid()

    # Step 3: Compute initial DEM
    est = x._estimate_initial_dem()

    # TEST 1: Correct length
    assert len(est) == len(x.logT)

    # TEST 2: All values should be exactly 0.0 ( Python implementation overrides with flat logDEM = 0)
    assert np.allclose(est, 0.0)

    # TEST 3: Internal storage _initial_log_dem should match
    assert np.allclose(x._initial_log_dem, est)

    # TEST 4: Returned DEM should be finite
    assert np.all(np.isfinite(est))



def test_prepare_spline_system():
    # Setup: 3 channels → n_spl = 2
    
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([1000.0, 2000.0, 1500.0])
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=8.0,
        logarithmic_temperature_step_size=0.1,
    )

    # Need logT, response matrix, initial DEM
    x.create_logT_grid()
    x._interpolate_responses_to_grid()
    x._estimate_initial_dem()   # sets _initial_log_dem

    # Prepare spline system
    x._prepare_spline_system()


    # TEST 1 — n_spl formula: n_channels=3 → n_spl=2
    assert x.n_spl == 2

    # TEST 2 — spline_logT: correct shape and increasing
    assert len(x.spline_logT) == 2
    assert np.all(np.diff(x.spline_logT) > 0)

    # TEST 3 — spline_log_dem should be zeros (flat initial DEM)
    assert len(x.spline_log_dem) == 2
    assert np.allclose(x.spline_log_dem, 0.0)

    # TEST 4 — pm_matrix has correct shape
    assert x.pm_matrix.shape == (3, len(x.logT))

    # pm_matrix should be >= 0
    assert np.all(x.pm_matrix >= 0)

    # TEST 5 — weights and abundances
    assert np.all(x.weights == 1.0)
    assert np.all(x.abundances == 1.0)

# def test_interpolate_responses_to_grid():
#     # -------------------------------
#     # Step 1: Setup a simple DEM case
#     # -------------------------------
#     filters = ["Al-poly", "Ti-poly"]
#     intensities = np.array([1000.0, 2000.0])

#     responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

#     x = XRTDEMIterative(
#         observed_channel=filters,
#         observed_intensities=intensities,
#         temperature_responses=responses,
#         minimum_bound_temperature=5.5,
#         maximum_bound_temperature=8.0,
#         logarithmic_temperature_step_size=0.1,
#     )

#     # -------------------------------
#     # Step 2: Create temperature grid
#     # -------------------------------
#     x.create_logT_grid()

#     # -------------------------------
#     # Step 3: Interpolate responses
#     # -------------------------------
#     x._interpolate_responses_to_grid()

#     # -------------------------------
#     # TEST 1: Correct response matrix shape
#     # -------------------------------
#     n_filters = len(filters)
#     n_T = len(x.logT)

#     assert x.response_matrix.shape == (n_filters, n_T)
#     # TEST 2: Response values must be non-negative
#     assert np.all(x.response_matrix >= 0)

#     # TEST 3: Responses at the boundaries are very small but not negative
#     assert np.all(x.response_matrix[:, 0] >= 0)
#     assert np.all(x.response_matrix[:, -1] >= 0)

#     # TEST 4: Values near boundaries should be small
#     assert np.all(x.response_matrix[:, 0] < 1e-27)
#     assert np.all(x.response_matrix[:, -1] < 1e-27)


#     # TEST 3: Interpolated values are finite

#     assert np.all(np.isfinite(x.response_matrix))

#     # TEST 4: Boundary values should be significantly smaller than peak response
#     for i in range(n_filters):
#         peak = np.max(x.response_matrix[i])
#         left = x.response_matrix[i, 0]
#         right = x.response_matrix[i, -1]

#         # Boundaries should be at least 100x smaller than peak
#         assert left < peak / 100
#         assert right < peak / 100

#         # -------------------------------
#         # TEST 5: Temperature ordering preserved
#         # Response_matrix row i should roughly follow original shape:
#         # no reversed ordering, no nan blocks.
#         # -------------------------------
#         for i in range(n_filters):
#             # Check no negative values
#             assert np.all(x.response_matrix[i] >= 0)

#             # Should have at least 2 non-zero values inside range
#             assert np.count_nonzero(x.response_matrix[i]) > 2

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

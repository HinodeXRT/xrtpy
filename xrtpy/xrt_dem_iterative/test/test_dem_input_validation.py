from pathlib import Path

import astropy.units as u
import numpy as np
import pytest
from lmfit import Parameters
from scipy.io import readsav
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


def test_validate_inputs_accepts_typical_three_filter_case():
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([500.0, 1800.0, 1820], dtype=float)
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    assert len(responses) == len(filters)

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=8.0,
        logarithmic_temperature_step_size=0.1,
        monte_carlo_runs=0,
    )

    x.validate_inputs()


def test_create_logT_grid_respects_bounds_and_step():
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([500.0, 1800.0, 1820], dtype=float)
    responses = generate_temperature_responses(filters, "2020-10-17T03:40:10")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.1,
    )

    x.create_logT_grid()

    assert (x.logT[0], x.logT[-1], x.dlogT) == pytest.approx((5.5, 7.5, 0.1))
    assert len(x.logT) == 21  # (7.5-5.5)/0.1 + 1 = 21
    assert np.isclose(x.dlnT, np.log(10) * 0.1)


def test_validate_inputs_accepts_defaults():
    filters = ["Be-thin", "Be-med"]
    i_obs = [100.0, 200.0]
    resp = generate_temperature_responses(filters, "2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp)
    dem.validate_inputs()  # Should NOT raise


def test_validate_inputs_rejects_mismatched_intensity_errors():
    filters = ["Be-thin", "Be-med"]
    i_obs = [100.0, 200.0]
    i_err = [20.0]  # Wrong length - should be two error/ uncertainties
    resp = generate_temperature_responses(filters, "2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp, intensity_errors=i_err)

    with pytest.raises(ValueError, match="intensity_errors must match"):
        dem.validate_inputs()


def test_create_logT_grid():

    filters = ["Al-poly/Ti-poly","C-poly","Be-thin"]
    intensities = np.array([150.0,670.0,392.0])
    responses = generate_temperature_responses(filters, "2018-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.1,
    )

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


def test_estimate_initial_dem_returns_flat_log_dem_zero():
    filters = ["Al-poly", "Ti-poly","Al-thick"]
    intensities = np.array([150.2, 230.0, 1321.1])
    responses = generate_temperature_responses(filters, "2012-10-27T23:02:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
    )

    # Step 2: Create temperature grid, response matrix, and initial DEM
    x.create_logT_grid()
    x._interpolate_responses_to_grid()
    est = x._estimate_initial_dem()

    # TEST 1: Correct length
    assert len(est) == len(x.logT)

    # TEST 2: All values should be exactly 0.0 ( Python implementation overrides with flat logDEM = 0)
    assert np.allclose(est, 0.0)

    # TEST 3: Internal storage _initial_log_dem should match
    assert np.allclose(x._initial_log_dem, est)

    # TEST 4: Returned DEM should be finite
    assert np.all(np.isfinite(est))


def test_prepare_spline_system_initializes_all_solver_state():
    """
    1. _prepare_spline_system runs without errors
    2. n_spl computed correctly
    3. spline_logT shape and monotonicity
    4. spline_log_dem has correct values
    5. pm_matrix has correct shape
    6. weights and abundances are all ones
    """
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([102.0, 202.9, 1500.1])
    responses = generate_temperature_responses(filters, "2012-10-27T12:30:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )


    x.create_logT_grid()
    x._interpolate_responses_to_grid()
    x._estimate_initial_dem()  # sets _initial_log_dem


    x._prepare_spline_system()

    # n_spl formula: n_channels=3 > n_spl=2
    assert x.n_spl == 2

    # spline knot positions
    assert len(x.spline_logT) == 2
    assert np.all(np.diff(x.spline_logT) > 0)

    # spline_log_dem should be zeros (flat initial DEM)
    assert len(x.spline_log_dem) == 2
    assert np.allclose(x.spline_log_dem, 0.0)

    # pm_matrix has correct shape
    assert x.pm_matrix.shape == (3, len(x.logT))

    # pm_matrix should be >= 0
    assert np.all(x.pm_matrix >= 0)

    # weights and abundances
    assert np.all(x.weights == 1.0)
    assert np.all(x.abundances == 1.0)


def test_residuals_matches_hand_computed_forward_model():
    """
    Create a fully synthetic DEM / response case so the forward model has a predictable value.
    This isolates and tests the math inside `_residuals`.
    """
    filters = ["Al-poly"]
    intensities = np.array([10.0])
    responses = generate_temperature_responses(filters, "2018-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=6.5,
        logarithmic_temperature_step_size=0.1,
    )

    x.create_logT_grid()
    n_T = len(x.logT)

    x.pm_matrix = np.full((1, n_T), 2.0)

    x.spline_logT = np.array([x.logT[0], x.logT[-1]])
    x.spline_log_dem = np.array([0.0, 0.0])
    x.n_spl = 2

    params = Parameters()
    params.add("knot_0", value=0.0, min=-20, max=0)
    params.add("knot_1", value=0.0, min=-20, max=0)

    x.intensities_scaled = np.array([10.0])
    x.sigma_scaled_intensity_errors = np.array([1.0])
    x.abundances = np.ones(1)
    x.weights = np.ones(1)

    residuals = x._residuals(params)

    expected = (2.0 * n_T - 10.0) / 1.0
    assert residuals.shape == (1,)
    assert residuals[0] == pytest.approx(expected)


def test_solve_single_dem_returns_zeros_when_all_intensities_zero():
    """
    If all observed intensities are zero, the solver must return:
        - DEM = all zeros or value of normalization_factor
        - modeled intensities = all zeros
        - chi sqr = 0
        - result = None
    This run will output a warning - expected due to no intensity_errors provided.
    """
    # filterwarnings = ignore:No intensity_errors provided

    filters = ["Al-poly", "Ti-poly","Al-mesh"]
    intensities = np.array([0.0, 0.0, 0])  # all zero → triggers nosolve
    responses = generate_temperature_responses(filters, "2022-11-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )

    x.create_logT_grid()
    x._interpolate_responses_to_grid()

    dem, modeled, chi2, result = x._solve_single_dem(
        observed_intensities_vals=intensities
    )

    # Assertions matching IDL behavior -DEM must be zero everywhere
    assert np.all(dem == 0.0)

    # modeled intensities must be all zero
    assert np.all(modeled == 0.0)

    # chi sqaure must be zero
    assert chi2 == 0.0

    # result object must be None (no lmfit minimization)
    assert result is None

    # Shape correctness
    assert dem.shape == x.logT.shape
    assert modeled.shape == intensities.shape

def test_monte_carlo_produces_non_identical_realizations(monkeypatch):
    """
    Monte Carlo DEM runs should produce DEMs that differ from the base DEM
    when observational noise causes perturbed intensities to change.

    This test verifies:
        -Monte Carlo output arrays have correct shapes
        -Each perturbed DEM differs from the base DEM
        -Perturbed intensities differ from the base intensities
        -No invalid values appear (no NaNs, no negatives)
    """
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([2300.0, 1500.0, 800.0], dtype=float)
    intensity_errors = 0.2 * intensities  # explicit sigma -> avoids warning and ensures perturbations
    responses = generate_temperature_responses(filters, "2012-10-27T08:23:54")
    n_runs = 5

    # Make Monte Carlo deterministic for CI/test stability
    _orig_default_rng = np.random.default_rng
    monkeypatch.setattr(np.random, "default_rng", lambda *a, **k: _orig_default_rng(0))

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        intensity_errors=intensity_errors,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.5,
        monte_carlo_runs=n_runs,
    )
    x.solve()

    # Shapes (contract)
    assert x.mc_dem.shape == (n_runs + 1, len(x.logT))
    assert x.mc_base_obs.shape == (n_runs + 1, len(filters))
    assert x.mc_mod_obs.shape == (n_runs + 1, len(filters))
    assert x.mc_chisq.shape == (n_runs + 1,)

    # MC should actually perturb observations and (usually) DEMs
    assert np.any(x.mc_base_obs[1:] != x.mc_base_obs[0])
    assert np.any(x.mc_dem[1:] != x.mc_dem[0])

    # No invalid values
    assert np.all(np.isfinite(x.mc_dem)) and np.all(x.mc_dem >= 0.0)
    assert np.all(np.isfinite(x.mc_mod_obs))
    assert np.all(np.isfinite(x.mc_chisq))

def test_reconstruct_dem_from_knots_matches_endpoints_and_is_finite():
    filters = ["Al-poly"]
    intensities = np.array([123.0], dtype=float)
    responses = generate_temperature_responses(filters, "2023-10-27T12:22:00")

    x = XRTDEMIterative(filters, intensities, responses)
    x.create_logT_grid()

    x.spline_logT = np.array([x.logT[0], x.logT[-1]])
    x.n_spl = 2
    x.spline_log_dem = np.array([-2.0, 1.0])

    params = Parameters()
    params.add("knot_0", value=-2.0)
    params.add("knot_1", value=1.0)

    dem = x._reconstruct_dem_from_knots(params)

    assert dem.shape == x.logT.shape
    assert np.isfinite(dem).all() and (dem >= 0.0).all()

    assert dem[0] == pytest.approx(1e-2, rel=1e-6)
    assert dem[-1] == pytest.approx(1e1, rel=1e-6)

    assert dem[0] < dem[len(dem) // 2] < dem[-1]


def test_solve_end_to_end_produces_finite_dem_and_mc_outputs():
    """
    Full DEM solving pipeline using real XRT filter responses.
    Verifies:
        -solve() runs without errors
        -base DEM exists, finite, and positive
        -modeled intensities computed
        -chi-square finite
        -Monte Carlo arrays created correctly (when N>0)
        -No NaNs, no negative DEM, no shape mismatches
    """

    filters = ["Al-poly", "Ti-poly", "Be-thin", "C-poly"]
    intensities = np.array([250.0, 1800.0, 1390.0, 453.0], dtype=float)
    responses = generate_temperature_responses(filters, "2021-02-17T23:59:30")
    N = 3

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=N,
    )

    x.solve()

    # Base DEM checks
    assert hasattr(x, "dem")
    assert x.dem.shape == x.logT.shape
    assert np.all(np.isfinite(x.dem))
    assert np.all(x.dem >= 0.0)

    # Peak must not be zero everywhere
    assert np.max(x.dem) > 0.0

    # Modeled intensities
    assert hasattr(x, "modeled_intensities")
    assert x.modeled_intensities.shape == intensities.shape
    assert np.all(np.isfinite(x.modeled_intensities))

    # Chi-square
    assert hasattr(x, "chisq")
    assert np.isfinite(x.chisq)

    # Monte Carlo arrays exist and are valid
    assert x.mc_dem.shape == (N + 1, len(x.logT))
    assert x.mc_base_obs.shape == (N + 1, len(filters))
    assert x.mc_mod_obs.shape == (N + 1, len(filters))
    assert x.mc_chisq.shape == (N + 1,)

    # MC fields finite
    assert np.all(np.isfinite(x.mc_dem))
    assert np.all(x.mc_dem >= 0.0)
    assert np.all(np.isfinite(x.mc_base_obs))
    assert np.all(np.isfinite(x.mc_mod_obs))
    assert np.all(np.isfinite(x.mc_chisq))

    # At least one MC DEM must differ from base DEM
    base_dem = x.mc_dem[0]
    different = [
        not np.allclose(base_dem, x.mc_dem[i], rtol=1e-5, atol=1e-8)
        for i in range(1, N + 1)
    ]
    assert any(different), "Monte Carlo DEMs identical to base DEM — noise not applied?"

    # Modes must differ for perturbed cases
    base_mod = x.mc_mod_obs[0]
    different_mod = [
        not np.allclose(base_mod, x.mc_mod_obs[i]) for i in range(1, N + 1)
    ]
    assert any(different_mod)

def test_warns_when_intensity_exceeds_physical_xrt_limit():
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([4200.0, 1800.0], dtype=float)
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(filters, intensities, responses)

    with pytest.warns(UserWarning, match=r"observed intensit"):
        x.validate_inputs()

def test_warns_when_intensity_is_negative():
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([-50.0, 1200.0], dtype=float)
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(filters, intensities, responses)

    with pytest.warns(UserWarning, match=r"negative"):
        x.validate_inputs()

# # ----------------------------------- TEST Against IDL ------------------------------------------


# def test_compare_with_idl_dem():
#     """
#     Compare xrtpy DEM solver output with reference IDL DEM stored in a .sav file.
#     The comparison is tolerant (allclose), since Python spline fitting and
#     MPFIT/LMFit differences mean we cannot expect bit-exact equality.
#     """

#     TEST_DIR = Path(__file__).parent
#     data_path = (
#         TEST_DIR / "IDL_DEM_testing_sav_files" / "obs_20090730_DEM_MC100_IDL_2026.sav"# "xrt_IDL_DEM_2012_10_27_MC100.sav"
#     )
#     data = readsav(data_path)
#     logT_idl = data["logt"]

#     # XRTpy - 2012-10-27 Data Set
#     # filters = ["Ti-poly", "Be-thin", "Al-poly", "C-poly"]
#     # intensities = np.array([311680.2, 135815.0, 2351258.9, 2352987.7])
#     # date = "2012-10-27T16:27:46"

#     filters = ["Al-mesh", "Ti-poly", "Al-poly", "Be-thin"]
#     intensities = [178.482 ,44.919,132.193,3.149]  # DN/s
#     observation_date="2009-07-30T00:38"

#     responses = generate_temperature_responses(filters, observation_date)

#     x = XRTDEMIterative(
#         observed_channel=filters,
#         observed_intensities=intensities,
#         temperature_responses=responses,
#         monte_carlo_runs=0,  # IDL comparison → no MC
#     )

#     # STEP 3 — Run Python DEM
#     x.solve()  # linear DEM
#     logT_python = x.logT

#     # STEP 4 — Compare logT grids - IDL usually uses exactly the same grid, but check with tolerance
#     assert np.allclose(
#         logT_python, logT_idl, atol=1e-6
#     ), "Temperature grids differ significantly between IDL and Python."
#############
# #log_dem_python = np.log10(np.maximum(dem_python, 1e-99))
# #dem_idl = data["dem"]
# #log_dem_idl = np.log10(np.maximum(dem_idl, 1e-99))
# # STEP 5 — Compare DEM values
# # Numerical tolerances:
# # - 0.3 dex tolerance → factor of ~2
# # - Acceptable because:
# #   • Spline fits differ slightly between MPFIT and LMFit
# #   • IDL interpolation behavior differs from CubicSpline
# #   • Floating-point roundoff differences
# tol_dex = 0.9

# diff = np.abs(log_dem_python - log_dem_idl)

# assert np.all(diff < tol_dex), (
#     "DEM shape diverges too far from IDL reference.\n"
#     f"Max difference {diff.max():.3f} dex (allowed {tol_dex})."
# )

# # Ensure the overall trends match: peak in same region
# peak_idl = np.argmax(log_dem_idl)
# peak_python = np.argmax(log_dem_python)
# assert abs(peak_idl - peak_python) <= 1, \
#     "DEM peak location differs more than 1 temperature bin."

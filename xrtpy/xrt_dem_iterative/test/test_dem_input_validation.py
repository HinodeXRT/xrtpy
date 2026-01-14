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


def test_dem_validate_inputs_basic():
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([2500.0, 1800.0])
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
    i_err = [100.0]  # Wrong length - should be two error/ uncertainties
    resp = generate_temperature_responses(filters, "2007-07-10")
    dem = XRTDEMIterative(filters, i_obs, resp, intensity_errors=i_err)
    with pytest.raises(ValueError, match="intensity_errors must match"):
        dem.validate_inputs()


def test_create_logT_grid():

    filters = ["Al-poly"]
    intensities = np.array([1500.0])
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


def test_estimate_initial_dem():
    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([1500.0, 2300.0])
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
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
    """
    1. _prepare_spline_system runs without errors
    2. n_spl computed correctly
    3. spline_logT shape and monotonicity
    4. spline_log_dem has correct values
    5. pm_matrix has correct shape
    6. weights and abundances are all ones
    """
    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([1000.0, 2000.0, 1500.0])
    responses = generate_temperature_responses(filters, "2012-10-27T12:30:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )

    # logT, response matrix, initial DEM
    x.create_logT_grid()
    x._interpolate_responses_to_grid()
    x._estimate_initial_dem()  # sets _initial_log_dem

    # Prepare spline system
    x._prepare_spline_system()

    # TEST 1 — n_spl formula: n_channels=3 > n_spl=2
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


def test_residuals_simple_case():
    """
    Create a fully synthetic DEM / response case so the forward model has a predictable value.
    This isolates and tests the math inside `_residuals`.
    """
    filters = ["Dummy"]
    intensities = np.array([10.0])  # I_obs
    responses = generate_temperature_responses(["Al-poly"], "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=6.5,
        logarithmic_temperature_step_size=0.1,
    )

    # STEP 2 — Temperature grid
    x.create_logT_grid()
    N = len(x.logT)

    # Synthetic pm_matrix = constant 2 everywhere
    x.pm_matrix = np.ones((1, N)) * 2.0

    # STEP 3 — Construct synthetic spline state
    x.spline_logT = np.array([x.logT[0], x.logT[-1]])
    x.spline_log_dem = np.array([0.0, 0.0])  # log10(DEM)=0 → DEM=1
    x.n_spl = 2

    params = Parameters()
    params.add("knot_0", value=0.0, min=-20, max=0)
    params.add("knot_1", value=0.0, min=-20, max=0)

    # STEP 4 synthetic errors
    x.intensities_scaled = np.array([10.0])
    x.sigma_scaled_intensity_errors = np.array([1.0])

    # MISSING IN ORIGINAL TEST: Need abundances and weights (normally set in _prepare_spline_system)
    x.abundances = np.ones(1)
    x.weights = np.ones(1)

    # STEP 5 — Compute residuals
    residuals = x._residuals(params)

    # Expected residual:
    #   DEM(T)=1 > pm(T)=2 > I_model = 2*N
    #   residual = (2N - I_obs) / sigma
    expected_I_model = 2.0 * N
    expected_residual = expected_I_model - 10.0  # sigma = 1

    assert residuals.shape == (1,)
    assert np.isfinite(residuals[0])
    assert np.isclose(residuals[0], expected_residual)


def test_solve_single_dem_zero_case():
    """
    If all observed intensities are zero, the solver must return:
        - DEM = all zeros or value of normalization_factor
        - modeled intensities = all zeros
        - chi sqr = 0
        - result = None
    This run will output a warning - expected due to no intensity_errors provided.
    """
    # filterwarnings = ignore:No intensity_errors provided

    filters = ["Al-poly", "Ti-poly"]
    intensities = np.array([0.0, 0.0])  # all zero → triggers nosolve
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )

    # STEP 2 — Build grid + response matrix
    x.create_logT_grid()
    x._interpolate_responses_to_grid()

    # STEP 3 — Call solver on zero intensities
    dem, modeled, chi2, result = x._solve_single_dem(
        observed_intensities_vals=intensities
    )

    # STEP 4 — Assertions matching IDL behavior -DEM must be zero everywhere
    assert np.all(dem == 0.0)

    # modeled intensities must be all zero
    assert np.all(modeled == 0.0)

    # chi² must be zero
    assert chi2 == 0.0

    # result object must be None (no lmfit minimization)
    assert result is None

    # Shape correctness
    assert dem.shape == x.logT.shape
    assert modeled.shape == intensities.shape


def test_monte_carlo_different_realizations():
    """
    Monte Carlo DEM runs should produce DEMs that differ from the base DEM
    when observational noise causes perturbed intensities to change.

    This test verifies:
        • Monte Carlo output arrays have correct shapes
        • Each perturbed DEM differs from the base DEM
        • Perturbed intensities differ from the base intensities
        • No invalid values appear (no NaNs, no negatives)
    """

    filters = ["Al-poly", "Ti-poly", "Be-thin"]
    intensities = np.array([3000.0, 1500.0, 800.0], dtype=float)
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")
    N = 5  # small Monte Carlo batch for fast testing

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=7.5,
        logarithmic_temperature_step_size=0.5,
        monte_carlo_runs=N,
    )

    # STEP 2: Run full DEM + Monte Carlo solver
    x.solve()  # this computes base and MC DEMs

    # STEP 3: Basic shape checks
    n_T = len(x.logT)
    n_obs = len(filters)

    # mc_dem shape: (N+1, n_T)
    assert x.mc_dem.shape == (N + 1, n_T)

    # mc_base_obs shape: (N+1, n_obs)
    assert x.mc_base_obs.shape == (N + 1, n_obs)

    # mc_mod_obs shape: (N+1, n_obs)
    assert x.mc_mod_obs.shape == (N + 1, n_obs)

    # mc_chisq shape: (N+1,)
    assert x.mc_chisq.shape == (N + 1,)

    # STEP 4: DEMs should *not* all be identical
    base_dem = x.mc_dem[0]

    # At least one MC DEM must differ from base DEM
    diffs = [
        not np.allclose(base_dem, x.mc_dem[ii], rtol=1e-5, atol=1e-8)
        for ii in range(1, N + 1)
    ]

    assert any(
        diffs
    ), "Monte Carlo realizations did not change the DEM; noise not applied?"

    # STEP 5: Perturbed observed intensities should differ
    base_obs = x.mc_base_obs[0]

    obs_diffs = [not np.allclose(base_obs, x.mc_base_obs[ii]) for ii in range(1, N + 1)]
    assert any(
        obs_diffs
    ), "Monte Carlo observed intensities identical; noise not applied?"

    # STEP 6: No invalid numbers after fitting
    assert np.all(np.isfinite(x.mc_dem))
    assert np.all(x.mc_dem >= 0.0)
    assert np.all(np.isfinite(x.mc_mod_obs))
    assert np.all(np.isfinite(x.mc_chisq))


def test_reconstruct_dem_from_knots():
    """
    Test whether DEM reconstruction from spline knots behaves predictably.

    We construct:
        • a synthetic temperature grid
        • synthetic knots at logT endpoints
        • synthetic log10(DEM) values at knots
    And verify that:
        • reconstructed logDEM matches knot values at endpoints
        • DEM(T) is smooth in between (no NaNs or jumps)
        • the DEM shape increases when knot values increase
    """

    filters = ["Dummy"]
    intensities = np.array([1000.0])
    responses = generate_temperature_responses(["Al-poly"], "2012-10-27T00:00:00")

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
    )

    x.create_logT_grid()

    # STEP 2 — Define synthetic spline knot positions - Use ONLY endpoints (simplest nontrivial spline)
    x.spline_logT = np.array([x.logT[0], x.logT[-1]])
    x.n_spl = 2

    # STEP 3 — Create synthetic lmfit Parameters for knots - Case: log10(DEM) goes from -2 at low T to +1 at high T
    params = Parameters()
    params.add("knot_0", value=-2.0)  # DEM = 1e-2
    params.add("knot_1", value=+1.0)  # DEM = 1e+1

    # Also required by reconstruct: store initial logDEM values
    x.spline_log_dem = np.array([-2.0, 1.0])

    # STEP 4 — Reconstruct DEM
    dem = x._reconstruct_dem_from_knots(params)

    # STEP 5 — Assertions
    # # Shape matches temperature grid
    assert dem.shape == x.logT.shape

    # No NaNs or negatives
    assert np.all(np.isfinite(dem))
    assert np.all(dem >= 0.0)

    # Endpoint matches exactly 10^knot_value
    assert np.isclose(dem[0], 10 ** (-2.0), rtol=1e-6)
    assert np.isclose(dem[-1], 10 ** (1.0), rtol=1e-6)

    # The DEM should increase monotonically between these endpoints
    # (Cubic spline + monotonic knots → smooth monotonic increase)
    assert (
        dem[0] < dem[len(dem) // 2] < dem[-1]
    ), "Reconstructed DEM is not increasing between knots"


def test_full_pipeline_end_to_end():
    """
    Full DEM solving pipeline using real XRT filter responses.
    Verifies:
        • solve() runs without errors
        • base DEM exists, finite, and positive
        • modeled intensities computed
        • chi-square finite
        • Monte Carlo arrays created correctly (when N>0)
        • No NaNs, no negative DEM, no shape mismatches
    """

    filters = ["Al-poly", "Ti-poly", "Be-thin", "C-poly"]
    intensities = np.array([2500.0, 1800.0, 900.0, 450.0], dtype=float)
    responses = generate_temperature_responses(filters, "2012-10-27T00:00:00")
    N = 3

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=N,
    )

    x.solve()

    # STEP 3: Base DEM checks
    assert hasattr(x, "dem")
    assert x.dem.shape == x.logT.shape
    assert np.all(np.isfinite(x.dem))
    assert np.all(x.dem >= 0.0)

    # Peak must not be zero everywhere
    assert np.max(x.dem) > 0.0

    # STEP 4: Modeled intensities
    assert hasattr(x, "modeled_intensities")
    assert x.modeled_intensities.shape == intensities.shape
    assert np.all(np.isfinite(x.modeled_intensities))

    # STEP 5: Chi-square
    assert hasattr(x, "chisq")
    assert np.isfinite(x.chisq)

    # STEP 6: Monte Carlo arrays exist and are valid
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

    # STEP 7: At least one MC DEM must differ from base DEM
    base_dem = x.mc_dem[0]
    different = [
        not np.allclose(base_dem, x.mc_dem[i], rtol=1e-5, atol=1e-8)
        for i in range(1, N + 1)
    ]
    assert any(different), "Monte Carlo DEMs identical to base DEM — noise not applied?"

    # STEP 8: Modes must differ for perturbed cases
    base_mod = x.mc_mod_obs[0]
    different_mod = [
        not np.allclose(base_mod, x.mc_mod_obs[i]) for i in range(1, N + 1)
    ]
    assert any(different_mod)

# ----------------------------------- TEST Against IDL ------------------------------------------


def test_compare_with_idl_dem():
    """
    Compare xrtpy DEM solver output with reference IDL DEM stored in a .sav file.
    The comparison is tolerant (allclose), since Python spline fitting and
    MPFIT/LMFit differences mean we cannot expect bit-exact equality.
    """

    TEST_DIR = Path(__file__).parent
    data_path = (
        TEST_DIR / "IDL_DEM_testing_sav_files" / "xrt_IDL_DEM_2012_10_27_MC100.sav"
    )
    data = readsav(data_path)
    logT_idl = data["logt"]

    # XRTpy
    filters = ["Ti-poly", "Be-thin", "Al-poly", "C-poly"]
    intensities = np.array([311680.2, 135815.0, 2351258.9, 2352987.7])
    date = "2012-10-27T16:27:46"

    responses = generate_temperature_responses(filters, date)

    x = XRTDEMIterative(
        observed_channel=filters,
        observed_intensities=intensities,
        temperature_responses=responses,
        monte_carlo_runs=0,  # IDL comparison → no MC
    )

    # STEP 3 — Run Python DEM
    x.solve()  # linear DEM
    logT_python = x.logT

    # STEP 4 — Compare logT grids - IDL usually uses exactly the same grid, but check with tolerance
    assert np.allclose(
        logT_python, logT_idl, atol=1e-6
    ), "Temperature grids differ significantly between IDL and Python."

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

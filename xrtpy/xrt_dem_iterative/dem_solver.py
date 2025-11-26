__all__ = [
    "XRTDEMIterative",
]

# import pdb; pdb.set_trace()

import warnings

import astropy.units as u
import numpy as np
from lmfit import Parameters, minimize
from scipy.interpolate import interp1d

from xrtpy.util.filters import validate_and_format_filters
from xrtpy.xrt_dem_iterative import dem_plotting


class XRTDEMIterative:
    """
    Differential Emission Measure (DEM) solver for Hinode/XRT observations.

    This class implements a Python version of the IDL routine
    `xrt_dem_iterative2.pro`, using spline-parameterized DEM curves and
    iterative least-squares fitting. It supports Monte Carlo error analysis
    and closely mirrors the logic of the original IDL algorithm.

    Parameters
    ----------
    observed_channel : str or list of str, required
        Names of the filters used in the observation (for example,
        "Al-mesh", "Be-thin"). Must correspond one-to-one with the
        temperature_responses argument.
    observed_intensities : array-like, required
        Observed intensities for each filter channel. Units are DN/s/pix.
    temperature_responses : list, required
        List of TemperatureResponseFundamental objects matching the filters.
        Units are DN s^-1 pix^-1 cm^5. These can be created using
        xrtpy.response.tools.generate_temperature_responses().
    intensity_errors : array-like, optional
        Uncertainties in the observed intensities. If None, a default model
        is used: max(0.03 * intensity, 2 DN/s/pix).
    minimum_bound_temperature : float,  optional
        Minimum value of the log10(T) grid. Default is 5.5.
    maximum_bound_temperature : float,  optional
        Maximum value of the log10(T) grid. Default is 8.0.
    logarithmic_temperature_step_size : float, optional
        Step size for the log10(T) grid. Default is 0.1.
    monte_carlo_runs : int, optional
        Number of Monte Carlo repetitions to perform. Default is 0 (disabled).
    max_iterations : int, optional
        Maximum number of function evaluations for lmfit. Default is 2000.
    normalization_factor : float, optional
        Internal scaling factor used during optimization. Default is 1e21.

    Notes
    -----
    - All lists (observed_channel, observed_intensities,
    temperature_responses) must be the same length.
    - The log10(T) range must lie inside the native temperature grid
    provided by all filter responses.
    - If intensity_errors is not provided, a default model is used to
    estimate uncertainties.

    SELFNOTEJOY
        Add web-link to IDL script.
    """

    def __init__(
        self,
        observed_channel,
        observed_intensities,
        temperature_responses,
        intensity_errors=None,
        minimum_bound_temperature=5.5,
        maximum_bound_temperature=8.0,
        logarithmic_temperature_step_size=0.1,
        monte_carlo_runs=0,
        max_iterations=2000,
        normalization_factor=1e21,
    ):
        """
        Notes
        -----
        - All input lists (`observed_channel`, `observed_intensities`, and `temperature_responses`)
        must be the same length. Each entry should correspond to one filter.

        - The temperature grid range (`minimum_bound_temperature`, `maximum_bound_temperature`) must lie entirely within the
        response temperature ranges for **all** filters provided.

        - If `intensity_errors` is not provided, a model-based error estimate is used:
        max(0.03 * observed_intensity, 2 (DN/s/pix)), as in the IDL original.

        - Default XRT filter names include:
        {'Al-mesh', 'Al-poly', 'C-poly', 'Ti-poly', 'Be-thin', 'Be-med', 'Al-med', 'Al-thick', 'Be-thick',
        'Al-poly/Al-mesh', 'Al-poly/Ti-poly', 'Al-poly/Al-thick', 'Al-poly/Be-thick'}
        """

        # Validate and store filter names
        if observed_channel is None or (
            hasattr(observed_channel, "__len__") and len(observed_channel) == 0
        ):
            raise ValueError("`observed_channel` is required and cannot be empty.")
        self.observed_channel = validate_and_format_filters(observed_channel)

        # Store intensity and error arrays
        if observed_intensities is None or len(observed_intensities) == 0:
            raise ValueError("`observed_intensities` is required and cannot be empty.")
        self._observed_intensities = np.asarray(observed_intensities, dtype=float)

        if not np.all(np.isfinite(self._observed_intensities)):
            raise ValueError("`observed_intensities` must be finite numbers.")

        # Errors
        # NEW (defer to validate_inputs)
        if intensity_errors is not None:
            self._intensity_errors = np.asarray(intensity_errors, dtype=float)
        else:
            self._intensity_errors = None

        # Store temperature grid parameters
        self._logarithmic_temperature_step_size = float(
            logarithmic_temperature_step_size
        )
        self._minimum_bound_temperature = float(minimum_bound_temperature)
        self._maximum_bound_temperature = float(maximum_bound_temperature)
        if not (self._minimum_bound_temperature < self._maximum_bound_temperature):
            raise ValueError(
                "minimum_bound_temperature must be < maximum_bound_temperature."
            )

        n_pts = (
            int(
                np.floor(
                    (self._maximum_bound_temperature - self._minimum_bound_temperature)
                    / logarithmic_temperature_step_size
                    + 1e-9
                )
            )
            + 1
        )
        if n_pts < 4:
            raise ValueError("Temperature grid must have at least 4 points.")

        # Validate Monte Carlo setting
        if isinstance(monte_carlo_runs, bool):
            raise TypeError(
                "monte_carlo_runs must be a non-negative whole number, not a boolean."
            )
        elif (
            isinstance(monte_carlo_runs, int | np.integer)
            or isinstance(monte_carlo_runs, float)
            and monte_carlo_runs.is_integer()
        ):
            self._monte_carlo_runs = int(monte_carlo_runs)
        else:
            raise ValueError(
                "monte_carlo_runs must be a non-negative whole number (e.g., 0, 1, 100). "
                "Decimal values are not allowed."
            )

        if self._monte_carlo_runs < 0:
            raise ValueError("monte_carlo_runs must be ≥ 0.")

        # Validate max_iterations
        if not isinstance(max_iterations, int | np.integer) or max_iterations <= 0:
            raise ValueError("max_iterations must be a positive integer.")

        self._max_iterations = int(max_iterations)

        # Check logarithmic_temperature_step_size is positive
        if self._logarithmic_temperature_step_size <= 0:
            raise ValueError(
                "logarithmic_temperature_step_size must be a positive scalar."
            )

        # Store temperature response objects
        self.responses = temperature_responses

        if temperature_responses is None or len(temperature_responses) == 0:
            raise ValueError("`temperature_responses` is required and cannot be empty.")

        # Validate that the temperature grid falls within the responses
        for r in self.responses:
            logT_grid = np.log10(r.temperature.to_value(u.K))
            if not (
                self._minimum_bound_temperature >= logT_grid.min()
                and self._maximum_bound_temperature <= logT_grid.max()
            ):
                raise ValueError(
                    f"The specified temperature range [{minimum_bound_temperature}, {maximum_bound_temperature}] is outside the bounds of one or more filter response grids.\n"
                    "Please ensure the temperature range fits within all responses.\n"
                    "Hint: Default response range is logT = 5.5 to 8.0. You can view each response's logT range via: [r.temperature for r in responses]"
                )

        # Check consistency between inputs
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

        try:
            value = float(normalization_factor)
        except (TypeError, ValueError) as err:
            raise ValueError(f"Invalid normalization_factor: {err}") from err

        if value <= 0:
            raise ValueError("normalization_factor must be a positive number.")

        self._normalization_factor = value
        # self._normalization_factor = value

        self._using_estimated_errors = (
            False  # track whether default error model has been used
        )

    def validate_inputs(self) -> None:
        """
        Validate user-provided inputs. Raises ValueError on any issue.

        Intentionally separate from __init__ so tests and users can construct
        the object first, then explicitly validate (matches test expectations).
        """
        # 1) observed_channel non-empty
        if self.observed_channel is None or len(self.observed_channel) == 0:
            raise ValueError("`observed_channel` is required and cannot be empty.")

        # 2) intensities: length & finite
        if self._observed_intensities is None or len(self._observed_intensities) == 0:
            raise ValueError("`observed_intensities` is required and cannot be empty.")
        if not np.all(np.isfinite(self._observed_intensities)):
            raise ValueError("`observed_intensities` must be finite numbers.")

        # 3) responses present
        if self.responses is None or len(self.responses) == 0:
            raise ValueError("`temperature_responses` is required and cannot be empty.")

        # 4) lengths consistent between filters/intensities/responses
        if not (
            len(self._observed_intensities)
            == len(self.responses)
            == len(self.observed_channel)
        ):
            raise ValueError(
                "Length mismatch: intensities, responses, and observed_channel must match."
            )

        # 5) temperature grid sanity
        if not (self._minimum_bound_temperature < self._maximum_bound_temperature):
            raise ValueError(
                "minimum_bound_temperature must be < maximum_bound_temperature."
            )
        if self._logarithmic_temperature_step_size <= 0:
            raise ValueError(
                "logarithmic_temperature_step_size must be a positive scalar."
            )
        n_pts = (
            int(
                np.floor(
                    (self._maximum_bound_temperature - self._minimum_bound_temperature)
                    / self._logarithmic_temperature_step_size
                    + 1e-9
                )
            )
            + 1
        )
        if n_pts < 4:
            raise ValueError("Temperature grid must have at least 4 points.")

        # 6) grid range inside every response
        for r in self.responses:
            # logT_grid = np√.log10(r.temperature.value)
            logT_grid = np.log10(r.temperature.to_value(u.K))
            if not (
                self._minimum_bound_temperature >= logT_grid.min()
                and self._maximum_bound_temperature <= logT_grid.max()
            ):
                raise ValueError(
                    f"The specified temperature range [{self._minimum_bound_temperature}, {self._maximum_bound_temperature}] "
                    "is outside the bounds of one or more filter response grids."
                )

        # 7) intensity_errors length & finiteness (only if provided)
        if self._intensity_errors is not None:
            if self._intensity_errors.shape != self._observed_intensities.shape:
                raise ValueError(
                    "Length of intensity_errors must match observed_intensities."
                )
            if not np.all(np.isfinite(self._intensity_errors)) or np.any(
                self._intensity_errors < 0
            ):
                raise ValueError("`intensity_errors` must be finite and >= 0.")

        # 8 warning
        if np.all(self._observed_intensities == 0):
            warnings.warn(
                "\n\n All observed intensities are zero. DEM solution will yield zero. "
                "Object created, but solving will return DEM=0. \n\n",
                stacklevel=2,
            )

        # success -> no return value
        return None

    def __repr__(self):
        return (
            f"<XRTDEMIterative(filters={self.filter_names}, "
            f"logT={self._minimum_bound_temperature:.2f}–{self._maximum_bound_temperature:.2f}, logarithmic_temperature_step_size={self._logarithmic_temperature_step_size:.3f})>"
        )

    # @property  #Removed if not used
    # def name(self) -> str:
    #     """
    #     The XRT filter channel name, standardized (e.g. "Al-mesh").
    #     """
    #     return self._name

    ##########################################################
    @property
    def observed_intensities(
        self,
    ) -> (
        u.Quantity
    ):  # Add method to account for known values not worth observed_intensities
        """
        Observed intensities with physical units.

        Returns
        -------
        `~astropy.units.Quantity`
            Intensities in DN/s/pix for each filter channel.
            Where "pix" means a one-arcsecond, full -resolution XRT pixel.
        """
        return self._observed_intensities * (u.DN / u.s)

    @property
    def filter_names(self):
        """
        List of filter names from the temperature responses.
        """
        return [r.filter_name for r in self.responses]

    @property
    def response_temperatures(self):
        """
        List of temperature grids (K) for each filter response.
        """
        return [r.temperature for r in self.responses]

    @property
    def response_values(self):
        """
        List of response values (DN cm^5 / pix / s) for each filter.
        """
        return [r.response for r in self.responses]

    @property
    def minimum_bound_temperature(self):
        """
        Lower bound of log10 temperature grid.
        """
        return self._minimum_bound_temperature

    @property
    def maximum_bound_temperature(self):
        """
        Upper bound of log10 temperature grid.
        """
        return self._maximum_bound_temperature

    @property
    def logarithmic_temperature_step_size(self):
        """
        Bin width of log10 temperature grid.
        """
        return self._logarithmic_temperature_step_size

    @property
    def min_observational_error(self):
        """
        Default - Minimum absolute observational intensity error applied to DN/s/pix when intensity error is not provided.
        """
        return 2 * (u.DN / u.s)

    @property
    def relative_error(self):
        """
        Relative error used to scale intensity if an error is not provided.
        Default is 0.03 (3%).
        """
        return 0.03

    @property
    def intensity_errors(self) -> u.Quantity:
        """
        Return the intensity uncertainty values.

        If the user supplied intensity_errors, those values are returned.
        Otherwise a default model is used:

            sigma = max(0.03 * intensity, 2 DN/s/pix)

        This behavior mirrors the default uncertainty logic of the IDL routine
        xrt_dem_iterative2.pro.

        `~astropy.units.Quantity`
        Intensity errors in DN/s for each filter.

        For details, see:
        https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro
        """
        if self._intensity_errors is not None:
            return self._intensity_errors * (u.DN / u.s)

        if not self._using_estimated_errors:
            warnings.warn(
                (
                    "\n\nNo intensity_errors provided. Using default model: "
                    "max(relative-error * observed_intensity, min_observational_error)\n"
                    f"=> relative_error = {self.relative_error}, "
                    f"min_observational_error = {self.min_observational_error.value} DN/s\n"
                    "See: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro\n\n"
                ),
                category=UserWarning,
                stacklevel=2,
            )

        self._using_estimated_errors = True

        # Fixed in units
        estimated = np.maximum(
            (self.relative_error * self._observed_intensities) * (u.DN / u.s),
            self.min_observational_error,
        )
        return estimated

    @property
    def monte_carlo_runs(self) -> int:
        """
        Return
        ------
        int
            Number of Monte Carlo runs to perform (0 = disabled).
        """
        return self._monte_carlo_runs

    @property
    def normalization_factor(self):
        """
        Scaling factor used during DEM optimization to stabilize the spline fit.
        Corresponds to `normalization_factor` in IDL (default 1e21).
        """
        return self._normalization_factor

    @property
    def max_iterations(self):
        """
        Maximum number of iterations allowed in the least-squares DEM solver
        (e.g., when using `lmfit.minimize`). Default is 2000.
        """
        return self._max_iterations

    def create_logT_grid(self):
        """
        Construct the regular log10 temperature grid for DEM calculations.

        This builds a regularly spaced grid in log10(temperature), then converts it
        to linear temperature for use in the DEM integral.

        Notes
        -----
        - IDL's `xrt_dem_iterative2.pro` describes this as the "regular logT grid".
        - Two forms of the temperature grid are stored:
            * self.logT : log10(T) values (dimensionless)
            * self.T    : linear temperatures (Kelvin, astropy.units.Quantity)
        - The grid is inclusive of both `minimum_bound_temperature` and `maximum_bound_temperature`, with step size `logarithmic_temperature_step_size`.

        Additional attributes created:
        - self.dlogT : float
            Step size in log10(T) (dimensionless).
        - self.dlnT : float
            Step size in natural log(T). Useful for IDL-style integrals of the form:
                F = int. DEM(T) * R(T) * T d(ln T)
        """
        # number of bins including endpoints - if default values are used - end value is 26
        self.n_bins = (
            int(
                round(
                    (self._maximum_bound_temperature - self._minimum_bound_temperature)
                    / self._logarithmic_temperature_step_size
                )
            )
            + 1
        )

        # This matches the IDL temperature grid exactly. self.logT & self.T.
        # inclusive logT grid (IDL-style regular grid)
        # Units = 'log K. Runs from minimum_bound_temperature to Max_T with bin-width = DT
        # SELFNOTEJOY- Do we need to add units - current holds no units- it's wokring correctly - Should this on the Test as well?- I don't think it- it's noted in IDL but used with units

        # np.linspace over np.arange - simple reproduces that reliably:endpoint included, exact number of bins, and no accumulating floating-point drift - Best match to IDL
        self.logT = np.linspace(
            self._minimum_bound_temperature,
            self._maximum_bound_temperature,
            self.n_bins,
        )

        # linear temperature grid in Kelvin
        self.T = (10.0**self.logT) * u.K

        self.dlogT = float(
            self._logarithmic_temperature_step_size
        )  # scalar spacing (dimensionless and natural-log equivalent)

        self.dlnT = (
            np.log(10.0) * self.dlogT
        )  # for IDL-style intergral DEM(T) * R(T) * T dlnT - IDL “regular logT grid”

    def _interpolate_responses_to_grid(self):
        """
        Interpolate all filter responses onto the solver's regular log10(T) grid.

        This constructs the response matrix used in the DEM forward model.
        Each filter's native temperature response (R(T)) is interpolated to the
        grid defined by self.logT. Extrapolated values outside the native
        response range are set to zero.

        Equivalent to the "Res_Mat" construction in the IDL routine xrt_dem_iterative2.pro.

        Notes
        -----
        - Response units (from XRTpy) are DN s^-1 pix^-1 cm^5.
        - Output matrix has shape (n_filters, n_temperatures).
        - Rows correspond to filters; columns correspond to temperature bins.

        -------
        IDL method of Interpolate emissivity.
        Interpolate all filter responses onto the common logT grid and build
        the response matrix.

        Equivalent to constructing `Res_Mat` in IDL's `xrt_dem_iterative2.pro`
        and in the DEM_Solver PDF documentation.

        Attributes Created
        ------------------
        interpolated_responses : list of ndarray
        _response_matrix : ndarray
            Stacked filter responses on the uniform logT grid.
        """
        if not hasattr(self, "logT"):
            raise AttributeError(
                "Temperature grid missing. Call create_logT_grid() first."
            )

        rows = []
        for T_orig, R_orig in zip(
            self.response_temperatures, self.response_values, strict=False
        ):  # Make sure that R_orig.value is indeed in DN/s/pix per cm^5
            logT_orig = np.log10(T_orig.to_value(u.K))

            response_vals = R_orig.value  # already in correct physical units for XRTpy #NOTEFORJOY- TRIPLE check this

            interp_func = interp1d(
                logT_orig,
                response_vals,
                kind="linear",
                bounds_error=False,
                fill_value=0.0,
                assume_sorted=True,
            )
            rows.append(interp_func(self.logT))

        self.interpolated_responses = rows
        self._response_matrix = np.vstack(rows).astype(float)

        # Store the physical unit for clarity
        # self._response_unit = u.DN / u.s / u.pix / (u.cm**5)
        self._response_unit = (u.DN / u.s / u.pix) * u.cm**5

        # Quick sanity check
        if self._response_matrix.shape != (len(self.responses), self.logT.size):
            raise RuntimeError("Interpolated response matrix has unexpected shape.")

    @property
    def response_matrix(self):
        """
        Response matrix (n_filters x n_temperatures) after interpolation.

        Units: DN s^-1 pix^-1 cm⁵ per emission measure.

        Equivalent to `Res_Mat` in IDL's `xrt_dem_iterative2.pro`.

        Raises
        ------
        AttributeError
            If `_interpolate_responses_to_grid()` has not been called yet.
        """
        if not hasattr(self, "_response_matrix"):
            raise AttributeError(
                "Response matrix not available. Call _interpolate_responses_to_grid() first."
            )
        return self._response_matrix

    def _prepare_scaled_observations(self):
        """
        Prepare the scaled observed intensities and uncertainties
        exactly as done in the IDL routine xrt_dem_iterative2.pro.

        IDL equivalent:
            input1.i_obs = input1.i_obs / solv_factor
            input1.i_err = input1.i_err / solv_factor
        """
        # Extract values as plain floats (DN/s/pix)
        intensities_scaled_raw = (
            self.observed_intensities.value
        )  # Might just remove this line and up in the normalization
        sigma_intensity_errors_raw = self.intensity_errors.to_value(
            u.DN / u.s
        )  # Might just remove this line and up in the normalization

        # Apply normalization
        self.intensities_scaled = intensities_scaled_raw / self.normalization_factor
        self.sigma_scaled_intensity_errors = (
            sigma_intensity_errors_raw / self.normalization_factor
        )

        # Store for solver
        self._scaled_prepared = True

    ############################ Everything line of code ABOVE is PREP for the DEM  #############################################

    # ****************************************************************************************************************************
    ############################ Everything line of code BELOW is FOR the DEM  ##################################################

    #############************************** Start of INITIAL ROUGH DEM ESTIMATE **************************##########################
    ################## An estimated EM shape based on simple intensity-over-response peaks, smoothed across T. #####################

    def _estimate_initial_dem(self, cutoff: float = 1.0 / np.e) -> np.ndarray:
        """
        Compute an initial DEM estimate closely following the structure of the IDL routine xrt_dem_iter_estim.

        The IDL code performs a rough DEM estimate by evaluating intensities
        relative to response peaks, but xrt_dem_iterative2 ultimately replaces
        that estimate with a flat log10(DEM) curve before calling the solver.

        This method repeats the peak-finding logic for diagnostic purposes, but
        the final DEM passed into the solver is always:

            log10(DEM(T)) = 0.0  for all temperature bins

        which corresponds to DEM(T) = 1 in arbitrary units. This reproduces the
        IDL initial condition exactly.

        Parameters
        ----------
        cutoff : float, optional
            Fraction of the peak response used to define the usable window
            around a channel's emissivity peak. Default is 1/e (approximately
            0.3679).

        Returns
        -------
        ndarray
            Initial log10(DEM) estimate on self.logT. This is always a flat
            array of zeros (IDL-equivalent behavior).
        """
        if not hasattr(self, "logT"):
            raise AttributeError(
                "Temperature grid missing. Call create_logT_grid() first."
            )
        if not hasattr(self, "_response_matrix"):
            raise AttributeError(
                "Response matrix missing. Call _interpolate_responses_to_grid() first."
            )

        # Optional: store the peak-based estimates for diagnostics only.
        # These are NOT used to set the initial DEM (IDL overwrites them
        # with a flat DEM before handing off to the solver).
        t_peaks = []
        log_dem_estimates = []

        # Loop over each filter/channel with non-zero intensity
        for T_orig, R_orig, I_obs in zip(
            self.response_temperatures,
            self.response_values,
            self._observed_intensities,
            strict=False,
        ):
            logT_orig = np.log10(T_orig.to_value(u.K))
            # Make sure the response is in DN s^-1 pix^-1 per EM (cm^-5)
            R_vals = R_orig.to_value((u.DN / u.s / u.pix) * u.cm**5)

            if I_obs <= 0 or np.all(R_vals <= 0):
                continue  # skip unusable channel

            # 1. Peak location (logT)
            max_idx = np.argmax(R_vals)
            peak_val = R_vals[max_idx]
            t_peak_raw = logT_orig[max_idx]

            # Round to nearest grid step in logT, similar to round_off(..., 0.1)
            step = self._logarithmic_temperature_step_size
            t_peak = np.round(t_peak_raw / step) * step

            # 2. Good window (where R > cutoff * peak)
            good = np.where(R_vals > peak_val * cutoff)[0]
            if good.size < 1:
                continue

            # 3. Compute denominator integral: sum(T * R * dlnT)
            T_good = 10.0 ** logT_orig[good]  # [K]
            R_good = R_vals[good]
            # Native spacing in log10(T)
            if logT_orig.size > 1:
                dlogT_native = np.mean(np.diff(logT_orig))
            else:
                # Degenerate case; fall back to solver grid spacing
                dlogT_native = step
            dlnT_native = np.log(10.0) * dlogT_native
            denom = np.sum(T_good * R_good * dlnT_native)

            if denom <= 0:
                continue

            # 4. DEM estimate at the peak (for diagnostics only)
            dem_peak = I_obs / denom  # [cm^-5 K^-1]
            if dem_peak <= 0 or not np.isfinite(dem_peak):
                continue

            log_dem_est = np.log10(dem_peak)
            t_peaks.append(t_peak)
            log_dem_estimates.append(log_dem_est)

        # Compact duplicate peak temperatures by averaging (diagnostic only)
        if t_peaks:
            uniq = {}
            for t, val in zip(t_peaks, log_dem_estimates, strict=False):
                uniq.setdefault(t, []).append(val)
            t_peaks_uniq = np.array(sorted(uniq.keys()))
            log_dem_uniq = np.array([np.mean(uniq[t]) for t in t_peaks_uniq])
            # Store raw estimated peaks for debugging/inspection if desired
            self._raw_estimated_dem_peaks = (t_peaks_uniq, log_dem_uniq)
        else:
            self._raw_estimated_dem_peaks = (np.array([]), np.array([]))

        # IDL BEHAVIOR: override with flat initial DEM
        # xrt_dem_iter_estim ultimately does: dem = 0.0*findgen(nt) + 1.0  ; Use flat dem for initial guess on a regular logT grid. We mirror that here exactly:

        # est_log_dem_on_grid = np.ones_like(self.logT, dtype=float) * 1.0 NOV20
        # est_log_dem_on_grid = np.ones_like(self.logT, dtype=float) * 0.0 #NOTEFORJOY
        est_log_dem_on_grid = np.zeros_like(self.logT)

        # Return the intial first guessed DEM

        # Store for later use by the solver
        self._initial_log_dem = est_log_dem_on_grid

        return est_log_dem_on_grid

    #############************************** End of INITIAL DEM ESTIMATE **************************##################################

    # -------------------------------------------------------------------------------------------------------------------------------

    def _prepare_spline_system(self):
        """
        Prepare the spline-based DEM parameterization.

        This mirrors the IDL routine mp_prep and sets up all arrays needed by
        the least-squares solver, including:

        - self.n_spl          : number of spline knots
        - self.spline_logT    : knot positions (evenly spaced in log10(T))
        - self.spline_log_dem : initial values of log10(DEM) at each knot
        - self.pm_matrix      : response matrix multiplied by T * d(ln T)
        - self.weights        : all ones (IDL uses a channel weighting mask)
        - self.abundances     : all ones

        pm_matrix corresponds to:
            pm[i, j] = R_i(T_j) * T_j * d(ln T)

        which appears in the forward model:
            I_model_i = sum_j DEM(T_j) * pm[i, j]
        """

        # Number of channels
        n_line = len(self._observed_intensities)

        # IDL: n_spl = min(n_line - 1, 7)  - Make this a keyword in the class so use can tune it?
        self.n_spl = min(max(n_line - 1, 1), 7)

        # Weights and abundances (IDL sets all =1)
        # Later, should I a use_line mask (IDL ignores lines with i_obs=0), but you can add that when you need it.
        self.weights = np.ones(n_line, dtype=float)
        self.abundances = np.ones(n_line, dtype=float)

        # pm_matrix = R(T) * T * dlnT     (IDL line: emis * 10^t * alog(10^dt))
        # units - DN/s/pix/cm^5 * K * dLnT * DEM == DN/s/PIX
        T_linear = self.T.to_value(u.K)

        self.pm_matrix = (self._response_matrix * T_linear * self.dlnT).astype(float)

        # Knot positions are evenly spaced in logT (IDL spl_t)
        self.spline_logT = np.linspace(self.logT.min(), self.logT.max(), self.n_spl)

        # Initial spline DEM values: sample from initial logDEM grid
        # (IDL spline(est_t, est_dem, spl_t))
        interp_init = interp1d(
            self.logT,
            self._initial_log_dem,  # IDL is flat logDEM = 1.0
            kind="linear",  # IDL uses a cubic spline later NOTEFORJOY NOV20
            bounds_error=False,
            fill_value="extrapolate",
        )

        self.spline_log_dem = interp_init(self.spline_logT)

    def _build_lmfit_parameters(self):
        """
        Build lmfit.Parameters object representing log10(DEM) at the spline knots.
        IDL limits each spline DEM parameter to [-20, 0].
        """

        if not hasattr(self, "spline_log_dem"):
            raise RuntimeError("Run _prepare_spline_system() first.")

        params = Parameters()

        for i in range(self.n_spl):
            params.add(
                f"knot_{i}",
                value=float(self.spline_log_dem[i]),
                min=-20.0,
                max=0.0,
                vary=True,
            )

        return params

    def _reconstruct_dem_from_knots(self, params):
        """
        Construct DEM(T) on self.logT using spline of log10(DEM) at knot positions.
        Uses a natural cubic spline interpolation in log10(DEM) space.
        """
        from scipy.interpolate import CubicSpline

        knot_vals = np.array([params[f"knot_{i}"].value for i in range(self.n_spl)])

        # Or used the code above but switch from linear to kind="cubic"
        cs = CubicSpline(self.spline_logT, knot_vals, bc_type="natural")
        log_dem = cs(self.logT)
        dem = 10.0**log_dem
        return dem

    def _residuals(self, params):
        """
        Compute residuals for use by the least-squares optimizer.

        Residuals are computed as:
            residual_i = (I_model_i - I_observed_i) / sigma_i

        where:
            I_model = (pm_matrix @ DEM) * abundances

        Returns
        -------
        ndarray
            Residuals for each filter channel.
        """

        # 1. DEM(T)
        dem = self._reconstruct_dem_from_knots(params)

        # 2. Modeled intensities (IDL: i_mod = (dem ## pm) * abunds)
        i_mod = (self.pm_matrix @ dem) * self.abundances

        # 3. Observed (scaled)
        y_scaled = self.intensities_scaled  # i_obs / solv_factor
        sigma_scaled = self.sigma_scaled_intensity_errors

        # 4. Residuals = (i_mod - y_obs) * weights / sigma
        residuals = (i_mod - y_scaled) * self.weights / sigma_scaled

        # chi^2 history, mostly for debugging
        chi2_val = np.sum(residuals**2)
        if not hasattr(self, "_iteration_chi2"):
            self._iteration_chi2 = []
        self._iteration_chi2.append(chi2_val)

        return residuals

    def _solve_single_dem(self, observed_intensities_vals: np.ndarray):
        """
        This method solves the DEM for one set of intensities only, without
        Monte Carlo perturbation.
        """

        nf = self._normalization_factor

        # 1. scaled obs/errors
        self.intensities_scaled = observed_intensities_vals / nf
        sigma_phys = self.intensity_errors.to_value(u.DN / u.s)
        self.sigma_scaled_intensity_errors = sigma_phys / nf

        # 2. trivial nosolve case
        if np.all(self.intensities_scaled == 0.0):
            dem_model = np.zeros_like(self.logT)
            dem_phys = dem_model * nf
            modeled_intensities_phys = np.zeros_like(observed_intensities_vals)
            return dem_phys, modeled_intensities_phys, 0.0, None

        # 3. initial guess (log10 DEM_model on grid)
        init_log_dem = self._estimate_initial_dem()  # flat ~ 1.0 in IDL
        self._initial_log_dem = init_log_dem

        # 4. spline system using that initial guess
        self._prepare_spline_system()
        params0 = self._build_lmfit_parameters()  # values = initial_log_dem at knots

        # 5. run minimizer
        result = minimize(
            self._residuals, params0, max_nfev=self._max_iterations
        )  # method='leastsq'

        # THIS is the critical part – use *result.params*, not params0 <<<
        dem_model = self._reconstruct_dem_from_knots(result.params)  # DEM_model(T)
        dem_phys = dem_model * nf

        i_mod_scaled = (self.pm_matrix @ dem_model) * self.abundances
        modeled_intensities_phys = i_mod_scaled * nf

        resid = self._residuals(result.params)
        chisq = float(np.sum(resid**2))

        return dem_phys, modeled_intensities_phys, chisq, result

    # -------------------------------------------------------------------------------------------------------------------------------

    def _run_monte_carlo(self):
        """
        Replicates IDL's Monte Carlo loop.
        Produces:
            - self.mc_dem           shape (n_T, N+1)
            - self.mc_base_obs      shape (n_obs, N+1)
            - self.mc_mod_obs       shape (n_obs, N+1)
            - self.mc_chisq         shape (N+1,)

        Note that N+1 rows means: row 0 = base case, rows 1..N = MC.
        """

        n_obs = len(self._observed_intensities)
        nT = len(self.logT)
        N = self._monte_carlo_runs

        # Prepare arrays
        mc_dem = np.zeros((nT, N + 1))
        mc_base = np.zeros((n_obs, N + 1))
        mc_mod = np.zeros((n_obs, N + 1))
        mc_chi = np.zeros(N + 1)

        # Base run first (IDL puts real data in column 0)
        dem = self.dem  # already scaled back by normalization
        mc_dem[:, 0] = dem
        mc_base[:, 0] = self._observed_intensities  # unscaled
        mc_mod[:, 0] = self.modeled_intensities  # unscaled
        mc_chi[0] = self.current_chi2

        # --- Run the MC loops ---
        rng = np.random.default_rng()  # like systime(1)

        for ii in range(1, N + 1):
            # Step 1: Perturbed intensities (scaled)
            perturbed = (
                self.intensities_scaled
                + rng.normal(size=n_obs) * self.sigma_scaled_intensity_errors
            )
            perturbed = np.clip(perturbed, 0, None)

            # Store unscaled in mc_base
            mc_base[:, ii] = perturbed * self.normalization_factor

            # If all zero → nosolve=True
            if np.all(perturbed == 0):
                mc_dem[:, ii] = 0.0
                mc_mod[:, ii] = 0.0
                mc_chi[ii] = 0.0
                continue

            # Step 2: assign perturbed intensities
            self.intensities_scaled = perturbed

            # Step 3: Rebuild spline system
            self._prepare_spline_system()

            # Step 4: Solve via lmfit
            params = self._build_lmfit_parameters()
            out = minimize(self._residuals, params, max_nfev=self.max_iterations)

            # Step 5: Reconstruct DEM
            dem = self._reconstruct_dem_from_knots(out.params)
            dem_scaled = dem * self.normalization_factor  # unscale
            mc_dem[:, ii] = dem_scaled

            # Step 6: Compute modeled intensities
            modeled = (self.pm_matrix @ dem) * self.abundances
            mc_mod[:, ii] = modeled * self.normalization_factor

            # Step 7: Compute chi-square
            resid = self._residuals(out.params)
            mc_chi[ii] = np.sum(resid**2)

        # store results
        self.mc_dem = mc_dem
        self.mc_base_obs = mc_base
        self.mc_mod_obs = mc_mod
        self.mc_chisq = mc_chi

    # -------------------------------------------------------------------------------------------------------------------------------

    def solve(self):
        """
        High-level DEM solver.

        Python analogue of IDL's xrt_dem_iterative2.pro:

        1. Validate inputs.
        2. Build the logT grid and interpolate temperature responses.
        3. Solve ONE base DEM using the original (unperturbed) intensities.
        4. If Monte Carlo is requested (monte_carlo_runs > 0), perform N
        perturbed solves by adding Gaussian noise to the base intensities.
        5. Store all outputs on the instance for later analysis/plotting.

        After calling solve(), the following attributes are defined:

        Base solution
        -------------
        logT : ndarray (n_T,)
            log10 temperature grid [K].
        dem : ndarray (n_T,)
            Best-fit DEM(T) in physical units [cm^-5 K^-1].
        chisq : float
            Chi-square of the base fit (sum of squared residuals).
        modeled_intensities : ndarray (n_channels,)
            Best-fit modeled intensities in [DN s^-1 pix^-1].
        _base_fit_result : lmfit.MinimizerResult
            Full lmfit result object for diagnostics.

        Monte Carlo products
        --------------------
        mc_dem : ndarray (N+1, n_T)
            DEM curves for base (row 0) and each Monte Carlo run (rows 1..N),
            in physical units [cm^-5 K^-1].
        mc_chisq : ndarray (N+1,)
            Chi-square values for base (index 0) and each MC run.
        mc_base_obs : ndarray (N+1, n_channels)
            Observed intensities [DN s^-1 pix^-1] for base + each MC run.
            Row 0 = original observation; rows 1..N = perturbed.
        mc_mod_obs : ndarray (N+1, n_channels)
            Modeled intensities [DN s^-1 pix^-1] corresponding to mc_dem.
        """

        # Validate inputs (IDL: argument checks near top)

        self.validate_inputs()

        # 1) Build logT grid and response matrix - IDL: regular logT grid + interpolated emissivities

        self.create_logT_grid()
        self._interpolate_responses_to_grid()

        # Base observed intensities in physical units [DN/s/pix]
        base_obs_phys = np.asarray(self._observed_intensities, dtype=float)

        # 2) Solve BASE DEM (unperturbed intensities) Corresponds to the first call to xrt_dem_iter_nowidget in IDL.

        dem_base, mod_base, chisq_base, base_result = self._solve_single_dem(
            observed_intensities_vals=base_obs_phys
        )

        # Store base solution
        self.logT_solution = self.logT.copy()  # alias
        self.dem = dem_base  # [cm^-5 K^-1]
        self.chisq = chisq_base  # chi-square
        self.modeled_intensities = mod_base  # [DN/s/pix]
        self._base_fit_result = base_result

        # 3) Allocate Monte Carlo arrays (IDL: base_obs, dem_out, chisq, mod_obs)

        n_T = self.logT.size
        n_ch = base_obs_phys.size
        N = self.monte_carlo_runs

        self.mc_dem = np.zeros((N + 1, n_T), dtype=float)
        self.mc_chisq = np.zeros((N + 1,), dtype=float)
        self.mc_base_obs = np.zeros((N + 1, n_ch), dtype=float)
        self.mc_mod_obs = np.zeros((N + 1, n_ch), dtype=float)

        # Row 0 = base solution (unperturbed)
        self.mc_dem[0, :] = dem_base
        self.mc_chisq[0] = chisq_base
        self.mc_base_obs[0, :] = base_obs_phys
        self.mc_mod_obs[0, :] = mod_base

        # 4) Monte Carlo loop

        if N > 0:
            rng = np.random.default_rng()  # like IDL's systime(1) seeding

            # Intensity errors in physical units [DN/s/pix]
            sigma_phys = self.intensity_errors.to_value(u.DN / u.s)

            for ii in range(1, N + 1):
                # Lightweight progress indicator
                if ii % max(1, N // 20) == 0:
                    print(f"  - Monte Carlo run {ii}/{N}")

                # 4a) Perturb intensities: I' = I + N(0, sigma), clipped at 0
                noise = rng.normal(loc=0.0, scale=sigma_phys, size=base_obs_phys.shape)
                obs_pert = base_obs_phys + noise
                obs_pert = np.maximum(obs_pert, 0.0)  # IDL: >0 to avoid negatives

                # 4b) Solve DEM for this perturbed realization
                dem_i, mod_i, chisq_i, _ = self._solve_single_dem(
                    observed_intensities_vals=obs_pert
                )

                # 4c) Store Monte Carlo results
                self.mc_dem[ii, :] = dem_i
                self.mc_chisq[ii] = chisq_i
                self.mc_base_obs[ii, :] = obs_pert
                self.mc_mod_obs[ii, :] = mod_i

        # 5) Return DEM for convenience
        return self.dem

    def summary(self):
        """
        Print a detailed, diagnostic summary of the DEM solver state.

        This provides:
            - Input observation details
            - Temperature grid configuration
            - Response matrix status
            - Spline system configuration
            - Base DEM fit results
            - Monte Carlo statistics (if available)
            - Available plotting helpers
        """

        print("\n" + "=" * 76)
        print("          XRTpy DEM Iterative — Solver Summary")
        print("=" * 76)

        # -----------------------------------------------------
        print("\nINPUT DATA")
        print("-" * 70)
        print(f" Filters: {self.filter_names}")
        print(
            f" Observed Intensities: {np.array(self._observed_intensities)}  DN/s/pix"
        )
        print(f" Number of channels: {len(self._observed_intensities)}")

        # Error model
        if self._intensity_errors is not None:
            print(" Intensity Errors: User-provided")
        else:
            print(" Intensity Errors: Auto-estimated (3% of I, min=2 DN/s)")

        print(f" Error values (DN/s): {self.intensity_errors.to_value('DN/s')}\n")

        # -----------------------------------------------------
        print("\nTEMPERATURE GRID")
        print("-" * 70)
        if hasattr(self, "logT"):
            print(f" logT range: {self.logT[0]:.2f}  to  {self.logT[-1]:.2f}")
            print(f" Number of temperature bins: {len(self.logT)}")
            print(f" logT (grid spacing): {self.dlogT:.3f}")
            print(f" lnT (natural log spacing): {self.dlnT:.3f}")
        else:
            print(" Grid has not been constructed (call solve()).")

        # -----------------------------------------------------
        print("\nRESPONSE MATRIX")
        print("-" * 70)
        if hasattr(self, "_response_matrix"):
            print(
                f" Matrix shape:    {self._response_matrix.shape}  (filters x T bins)"
            )
            print(f" Response units: {self._response_unit}")
        else:
            print(" Response matrix not constructed.")

        # -----------------------------------------------------
        print("\nSOLVER CONFIGURATION")
        print("-" * 70)
        print(f" Normalization factor:     {self.normalization_factor:.2e}")
        print(f" Max iterations:           {self.max_iterations}")
        print(f" Monte Carlo runs:         {self.monte_carlo_runs}")
        if hasattr(self, "n_spl"):
            print(f" Number of spline knots:   {self.n_spl}")
            print(f" Knot positions (logT):    {getattr(self, 'spline_logT', 'N/A')}")
        else:
            print(" Spline system not prepared yet.")
        # -----------------------------------------------------
        print("\nINITIAL DEM GUESS")
        print("-" * 70)
        if hasattr(self, "_initial_log_dem"):
            print(" Initial DEM assumption:   flat log10(DEM) (IDL-style)")
            print(f" First 5 bins (log10):     {self._initial_log_dem[:5]}")
        else:
            print(" Initial DEM has not been estimated.")

        # -----------------------------------------------------

        print("\nBASE DEM SOLUTION")
        print("-" * 70)
        if hasattr(self, "dem"):
            print(f" DEM shape:               {self.dem.shape}")
            print(f" First 5 DEM bins:        {self.dem[:5]}")
            print(f" log10(DEM) first 5:      {np.log10(self.dem[:5] + 1e-99)}")
            print(f" Chi-square:              {self.chisq:.4e}")
            print(f" Modeled intensities:     {self.modeled_intensities}")
        else:
            print(" No DEM solution computed yet (call solve()).")

        # -----------------------------------------------------

        print("\nMONTE CARLO ENSEMBLE")
        print("-" * 70)
        if hasattr(self, "mc_dem"):
            N = self.mc_dem.shape[0] - 1
            print(f" MC realizations:         {N}")
            if N > 0:
                median = np.median(self.mc_dem[1:], axis=0)
                p16, p84 = np.percentile(self.mc_dem[1:], [16, 84], axis=0)
                print(" MC DEM statistics (first T-bin):")
                print(f"   Median (first 5):      {median[:5]}")
                print(
                    f"   1x bounds (log10):     "
                    f"{np.log10(p16[0] + 1e-99):.2f} – {np.log10(p84[0] + 1e-99):.2f}"
                )
            else:
                print(" MC array allocated but N=0 (no Monte Carlo).")
        else:
            print(" No Monte Carlo results available.")

        # -----------------------------------------------------
        print("\nPLOTTING HELPERS")
        print("-" * 76)
        print(" • plot_dem()         – Base DEM only")
        print(" • plot_dem_mc()      – Base DEM + MC ensemble")
        print("\n" + "=" * 76 + "\n")


XRTDEMIterative.plot_dem = dem_plotting.plot_dem
XRTDEMIterative.plot_dem_mc = dem_plotting.plot_dem_mc

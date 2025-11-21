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
    Estimate the differential emission measure (DEM) from Hinode/XRT data
    using the iterative spline-based method.

    Parameters
    ----------
    observed_channel : str or list of str (required)
        Filter names used in the observation (e.g., 'Al-mesh', 'Be-thin').
        Must match the provided temperature responses.
    observed_intensities : array-like (required)
        Observed intensities for each channel.
        Units = DN/s/pix.
    temperature_responses : list (required)
        List of `TemperatureResponseFundamental` objects matching the filters.
        Units = DN s^-1 pix^-1 EM^-1.
        Can be generated using `xrtpy.response.tools.generate_temperature_responses`
        for one or more filters. See: https://xrtpy.readthedocs.io/en/latest/getting_started.html
    intensity_errors : array-like, optional
        Intensity uncertainties. If None, will use a model-based estimate.
    minimum_bound_temperature : float
        Minimum log10 temperature (default: 5.5).
    maximum_bound_temperature: float
        Maximum log10 temperature (default: 8.0).
    logarithmic_temperature_step_size : float
        Step size in log10 temperature space (default: 0.1).
    monte_carlo_runs : int, optional
        Number of Monte Carlo runs to perform (default: 0, disabled).
        Each run perturbs `observed_intensities` using `intensity_errors`
        as Gaussian sigma and re-solves the DEM.
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
            raise ValueError(
                "monte_carlo_runs must be a non-negative whole number, not a boolean."
            )
        elif (
            isinstance(monte_carlo_runs, (int, np.integer))
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
        if not isinstance(max_iterations, (int, np.integer)) or max_iterations <= 0:
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

        # I am commenting this out because it is redundant since I am defining it below again. I wanna be consistent in using the same logT below.- Remove after testing
        # self.logT = np.arange(
        #     self._minimum_bound_temperature,
        #     self._maximum_bound_temperature
        #     + self._logarithmic_temperature_step_size / 2,
        #     self._logarithmic_temperature_step_size,
        # )

        try:
            self._normalization_factor = float(normalization_factor)
            if self._normalization_factor <= 0:
                raise ValueError("normalization_factor must be a positive number.")
        except Exception as e:
            raise ValueError(f"Invalid normalization_factor: {e}")

        self._using_estimated_errors = (
            False  # track whether default error model has been used
        )

    #### TEST GIT CI TEST #####

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
                "Object created, but solving will return DEM=0. \n\n"
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
        Returns the intensity uncertainties, either user-provided or model-based.

        If not provided, errors are estimated using:
            max(0.03 * observed_intensity, 2 DN/s/pix)

        For details, see:
        https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro

        Returns
        -------
        `~astropy.units.Quantity`
            Intensity errors in DN/s for each filter.
        """
        if self._intensity_errors is not None:
            return self._intensity_errors * (u.DN / u.s)

        # if self._using_estimated_errors:
        #     warnings.warn(
        #         "No intensity_errors provided. Using default model: "
        #         f"max(relative-error * observed_intensity, min_observational_error)\n"
        #         f"=> relative_error = {self.relative_error} =, min_observational_error = {self.min_observational_error.value} DN/s\n"
        #         "See: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro",
        #         UserWarning,
        #     )
        # self._using_estimated_errors = True  # suppress future warnings

        if not self._using_estimated_errors:
            warnings.warn(
                "No intensity_errors provided. Using default model: "
                f"max(relative-error * observed_intensity, min_observational_error)\n"
                f"=> relative_error = {self.relative_error} =, min_observational_error = {self.min_observational_error.value} DN/s\n"
                "See: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro",
                UserWarning,
            )
        self._using_estimated_errors = True

        # NOTETOJOYWe can remove if no issues later
        # #No units - added in the return
        # estimated = np.maximum(
        #     self.relative_error * self._observed_intensities ,
        #     self.min_observational_error.value,
        # )
        # return estimated * (u.DN / u.s)

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
        IDL method of Interpolate emissivity.
        Interpolate all filter responses onto the common logT grid and build
        the response matrix.

        Equivalent to constructing `Res_Mat` in IDL's `xrt_dem_iterative2.pro`
        and in the DEM_Solver PDF documentation.

        Notes
        -----
        - Each filter's response is interpolated to `self.logT` (regular log10 grid).
        - Extrapolation beyond the native response grid is set to 0.0.
        - Units: DN s^-1 pix^-1 cm^5 (per emission measure).
        - Shape of `_response_matrix`: (n_filters, n_temperatures)
        Rows = filters, Columns = temperature bins.

        Attributes Created
        ------------------
        interpolated_responses : list of ndarray
            Interpolated response arrays for each filter.
        _response_matrix : ndarray
            Final stacked matrix (n_filters x n_temperatures).
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
            # response_vals = R_orig.to_value((u.cm**5 * u.DN) / (u.pix * u.s))
            # response_vals = R_orig.to_value(u.DN / u.s / u.pix / (u.cm**5))
            # response_vals = R_orig.to_value((u.DN / u.s / u.pix) * u.cm**5) Comment on Nov14
            # response_vals = R_orig.to_value(u.DN / u.s / u.pix / u.cm**5)
            response_vals = (
                R_orig.value
            )  # already in correct physical units for XRTpy #NOTEFORJOY- TRIPLE check this

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

    ################################################################################################################################
    ################################################################################################################################
    #################################### structure with all fields the DEM solver expects  ##########################################
    # 1 Temperature - self.logT , self.T
    # 2 Response - note - interpolated onto your logT grid. - self.interpolated_responses, self._response_matrix
    # 3 # of bins  - n_bins
    # 4 i_obs – self._observed_intensites - measured DN/s/pixel scaled by solv_factor
    # self.observed_intensities
    # 5 uncertainty on the intensity - Also scaled by solv_factor.  - self.intensity_errors self.normalization_factor
    # 6 units?

    ################################################################################################################################
    ################################################################################################################################
    ################################################################################################################################

    # ****************************************************************************************************************************
    ############################ Everything line of code BELOW is FOR the DEM  ##################################################
    # ****************************************************************************************************************************

    # -------------------------------------------------------------------------------------------------------------------------------

    #############************************** Start of INITIAL ROUGH DEM ESTIMATE **************************##########################
    ################## An estimated EM shape based on simple intensity-over-response peaks, smoothed across T. #####################

    def _estimate_initial_dem(self, cutoff: float = 1.0 / np.e) -> np.ndarray:
        """
        Construct an initial DEM guess, mirroring IDL's xrt_dem_iter_estim.

        This method follows the *structure* of the IDL routine:
        - Identify channels with non-zero observed intensity.
        - For each such channel, find the peak of its emissivity/response.
        - Integrate the response around the peak to estimate a DEM value.
        - Combine/compact duplicate peak temperatures.

        HOWEVER, to exactly match the behavior of IDL's xrt_dem_iter_estim
        as used by xrt_dem_iter_nowidget, the final initial guess returned
        to the solver is a *flat* log10(DEM) curve:

            log10(DEM(T)) = 1.0  for all T on the solver grid.

        The detailed peak-based DEM estimates are kept only for optional
        diagnostics; they do not affect the initial DEM passed into the
        spline/least-squares solver (this is exactly what the IDL code does).

        Parameters
        ----------
        cutoff : float, optional
            Fraction of the peak response used to define the "good" window
            around each filter's peak. Default is 1/e (≈ 0.3679), as in IDL.

        Returns
        -------
        est_log_dem_on_grid : ndarray
            Array of shape (n_temperatures,) giving the initial guess for
            log10(DEM) on `self.logT`. For strict IDL-compatibility, this
            is identically 1.0 everywhere.
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
        # xrt_dem_iter_estim ultimately does:
        #   dem = 0.0*findgen(nt) + 1.0  ; Use flat dem for initial guess
        # on a regular logT grid. We mirror that here exactly:

        # est_log_dem_on_grid = np.ones_like(self.logT, dtype=float) * 1.0 NOV20

        # est_log_dem_on_grid = np.ones_like(self.logT, dtype=float) * 0.0 #NOTEFORJOY
        # or
        est_log_dem_on_grid = np.zeros_like(self.logT)

        # Return the intial first guessed DEM

        # Store for later use by the solver
        self._initial_log_dem = est_log_dem_on_grid

        return est_log_dem_on_grid

    #############************************** End of INITIAL DEM ESTIMATE **************************##################################

    # -------------------------------------------------------------------------------------------------------------------------------

    #############************************** Start of  **************************##########################

    def _prepare_spline_system(self):
        """
        Pythonic, IDL version of mp_prep.
        Prepares:s
            - self.n_spl           (number of spline knots)
            - self.spline_logT     (knot positions)
            - self.spline_log_dem  (initial spline logDEM values)
            - self.pm_matrix       (R(T) * T * dlnT)
            - self.weights         (all ones)
            - self.abundances      (all ones)
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
            kind="linear",  # IDL uses a cubic spline later NOTEFORJOY NOV24
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
        """
        from scipy.interpolate import CubicSpline

        knot_vals = np.array([params[f"knot_{i}"].value for i in range(self.n_spl)])

        # interp_spline = interp1d(
        #     self.spline_logT,
        #     knot_vals,
        #     kind="linear", #IDL uses cubic spline interpolation NOTEFORJOY NOV20
        #     bounds_error=False,
        #     fill_value="extrapolate",
        # )

        # log_dem = interp_spline(self.logT)  # log10(DEM)
        # dem = 10.0**log_dem  # DEM in linear cm^-5 K^-1

        # Or used the code above but switch from linear to kind="cubic"
        cs = CubicSpline(self.spline_logT, knot_vals, bc_type="natural")
        log_dem = cs(self.logT)
        dem = 10.0**log_dem
        return dem

    def _residuals(self, params):
        """
        IDL equivalent of mpdemfunct.
        Computes residuals:
            ((DEM ## pm) - i_obs_scaled) / i_err_scaled
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

        # Store χ² if desired
        chi2_val = np.sum(residuals**2)
        if not hasattr(self, "_iteration_chi2"):
            self._iteration_chi2 = []
        self._iteration_chi2.append(chi2_val)

        return residuals

    def _solve_single_dem(self, observed_intensities_vals: np.ndarray):
        """
        Solve the DEM once for a given set of observed intensities (base or MC-perturbed).

        Parameters
        ----------
        observed_intensities_vals : ndarray
            1D array of observed intensities in DN/s/pix (no units attached).

        Returns
        -------
        dem_phys : ndarray
            DEM(T) in physical units [cm^-5 K^-1] on self.logT.
        modeled_intensities_phys : ndarray
            Modeled intensities in DN/s/pix for each channel.
        chisq : float
            Sum of squared residuals for this run.
        fit_result : lmfit.MinimizerResult
            Full lmfit result object (for diagnostics).
        """
        # 1. Set scaled observations for this run (IDL: i_obs = i_obs/solv_factor)
        nf = self._normalization_factor
        self.intensities_scaled = observed_intensities_vals / nf
        # Errors are the same for all runs; scale once
        sigma_phys = self.intensity_errors.to_value(u.DN / u.s)
        self.sigma_scaled_intensity_errors = sigma_phys / nf

        # 2. If all intensities are zero → nosolve, DEM = 0
        if np.all(self.intensities_scaled == 0.0):
            dem_scaled = np.zeros_like(self.logT, dtype=float)
            dem_phys = dem_scaled * nf
            modeled_intensities_phys = np.zeros_like(observed_intensities_vals)
            chisq = 0.0
            fit_result = None
            return dem_phys, modeled_intensities_phys, chisq, fit_result

        # 3. Initial DEM & spline system (IDL: xrt_dem_iter_estim + mp_prep)
        self._estimate_initial_dem()
        self._prepare_spline_system()
        params = self._build_lmfit_parameters()

        # 4. Run the least-squares solver (IDL: xrt_dem_iter_solver + mpfit)
        result = minimize(self._residuals, params, max_nfev=self._max_iterations)

        # 5. Reconstruct DEM in *scaled* units, then convert to physical
        dem_scaled = self._reconstruct_dem_from_knots(result.params)  # cm^-5 K^-1 / nf
        dem_phys = dem_scaled * nf  # undo normalization, like IDL

        # 6. Modeled intensities (IDL: i_mod = dem ## pm * abunds)
        i_mod_scaled = (self.pm_matrix @ dem_scaled) * self.abundances
        modeled_intensities_phys = i_mod_scaled * nf  # back to DN/s/pix

        # 7. χ² from residuals
        resid = self._residuals(result.params)
        chisq = float(np.sum(resid**2))

        return dem_phys, modeled_intensities_phys, chisq, result

    ################################################################################################################################

    # -------------------------------------------------------------------------------------------------------------------------------

    #############************************** Start of  error bars / Monte Carlo  **************************##########################

    def _run_monte_carlo(self, result_params):
        """
        Replicates IDL's Monte Carlo loop.
        Produces:
            - self.mc_dem           shape (n_T, N+1)
            - self.mc_base_obs      shape (n_obs, N+1)
            - self.mc_mod_obs       shape (n_obs, N+1)
            - self.mc_chisq         shape (N+1,)
        """

        n_obs = len(self._observed_intensities)
        nT = len(self.logT)
        N = self._monte_carlo_runs

        # Prepare arrays
        mc_dem = np.zeros((nT, N + 1))
        mc_base = np.zeros((n_obs, N + 1))
        mc_mod = np.zeros((n_obs, N + 1))
        mc_chi = np.zeros(N + 1)

        # --- Base run first (IDL puts real data in column 0) ---
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

    #############************************** Start of DEM SOLVER  **************************##########################

    def solve(self):
        """
        High-level DEM solver.

        Replicates IDL’s xrt_dem_iterative2.pro behavior:
        1. Validate inputs
        2. Prepare logT grid and interpolate responses
        3. Solve ONE base DEM using original intensities
        4. If Monte Carlo requested, perform N perturbed solves
        5. Store all arrays cleanly for plotting and analysis

        After calling solve(), the following attributes exist:

        Base solution:
            - self.logT        (temperature grid)
            - self.dem         (DEM(T) in cm^-5 K^-1)
            - self.chisq       (chi-square of base fit)
            - self.modeled_intensities

        Monte Carlo products (N = monte_carlo_runs, N>=0):
            - self.mc_dem        shape = (N+1, n_T)
            - self.mc_chisq      shape = (N+1,)
            - self.mc_base_obs   shape = (N+1, n_filters)
            - self.mc_mod_obs    shape = (N+1, n_filters)

        Column 0 always holds the BASE solution (unperturbed).
        Columns 1..N hold Monte Carlo solutions.
        """

        # 0) Validate inputs -------------------------------------------------------
        self.validate_inputs()

        # 1) Build temperature grid and response matrix --------------------------
        self.create_logT_grid()
        self._interpolate_responses_to_grid()

        # Base observed intensities (physical DN/s/pix)
        base_obs_phys = self._observed_intensities.astype(float)

        # 2) Solve base DEM -------------------------------------------------------
        dem_base, mod_base, chisq_base, base_result = self._solve_single_dem(
            observed_intensities_vals=base_obs_phys
        )

        # Store base solution
        self.logT_solution = self.logT.copy()
        self.dem = dem_base
        self.chisq = chisq_base
        self.modeled_intensities = mod_base

        # 3) Allocate Monte Carlo arrays ------------------------------------------
        n_T = self.logT.size
        n_ch = base_obs_phys.size
        N = self.monte_carlo_runs

        self.mc_dem = np.zeros((N + 1, n_T), dtype=float)
        self.mc_chisq = np.zeros((N + 1,), dtype=float)
        self.mc_base_obs = np.zeros((N + 1, n_ch), dtype=float)
        self.mc_mod_obs = np.zeros((N + 1, n_ch), dtype=float)

        # Column 0 = base solution
        self.mc_dem[0, :] = dem_base
        self.mc_chisq[0] = chisq_base
        self.mc_base_obs[0, :] = base_obs_phys
        self.mc_mod_obs[0, :] = mod_base

        # 4) Monte Carlo Loop -----------------------------------------------------
        if N > 0:
            rng = np.random.default_rng()
            sigma_phys = self.intensity_errors.to_value(u.DN / u.s)

            for ii in range(1, N + 1):

                if ii % max(1, N // 20) == 0:  # print 5% updates
                    print(f"  - MC run {ii}/{N}")

                # Perturb intensities (IDL: i_obs + randn * i_err)
                noise = rng.normal(loc=0.0, scale=sigma_phys, size=base_obs_phys.shape)
                obs_pert = base_obs_phys + noise
                obs_pert = np.maximum(obs_pert, 0.0)  # IDL clips at zero

                # Solve DEM for perturbed intensities
                dem_i, mod_i, chisq_i, _ = self._solve_single_dem(obs_pert)

                # Store results
                self.mc_dem[ii, :] = dem_i
                self.mc_chisq[ii] = chisq_i
                self.mc_base_obs[ii, :] = obs_pert
                self.mc_mod_obs[ii, :] = mod_i

        # 5) Return DEM for convenience ------------------------------------------
        return self.dem


# Attach plotting functions from plotting.py to the class
XRTDEMIterative.plot_dem_results = dem_plotting.plot_dem_results
# XRTDEMIterative.plot_dem_uncertainty = dem_plotting.plot_dem_uncertainty
XRTDEMIterative.plot_idl_style = dem_plotting.plot_idl_style
XRTDEMIterative.plot_fit_residuals = dem_plotting.plot_fit_residuals
XRTDEMIterative.plot_dem_with_median_bins = dem_plotting.plot_dem_with_median_bins
XRTDEMIterative.plot_iteration_stats = dem_plotting.plot_iteration_stats
XRTDEMIterative.plot_dem = dem_plotting.plot_dem
XRTDEMIterative.plot_dem_uncertainty = dem_plotting.plot_dem_uncertainty
XRTDEMIterative.plot_mc_fan = dem_plotting.plot_mc_fan
XRTDEMIterative.plot_mc_hist_at_temperature = dem_plotting.plot_mc_hist_at_temperature
XRTDEMIterative.plot_mc_chisq = dem_plotting.plot_mc_chisq
XRTDEMIterative.plot_modeled_vs_observed = dem_plotting.plot_modeled_vs_observed

XRTDEMIterative.plot_dem_with_mc_vertical_bars = (
    dem_plotting.plot_dem_with_mc_vertical_bars
)
XRTDEMIterative.plot_dem_mc_vertical_bars = dem_plotting.plot_dem_mc_vertical_bars


################################################################################################################################
################################################################################################################################
#############************************** END of  DEM SOLVER **************************##########################
################################################################################################################################
################################################################################################################################


# def fit_with_multiple_methods(
#     self, methods=("leastsq", "least_squares", "nelder"), n_knots: int = 6, **kwargs
# ):
#     """
#     Try multiple lmfit minimization methods and pick the best χ².

#     Parameters
#     ----------
#     methods : tuple of str, optional
#         Minimization methods to test. Default = ("leastsq", "least_squares", "nelder").
#     n_knots : int, optional
#         Number of spline knots for DEM fit. Default = 6.
#     **kwargs : dict
#         Extra arguments passed to `lmfit.minimize`.

#     Returns
#     -------
#     best_result : lmfit.MinimizerResult
#         Result from the method with lowest chi².
#     """

#     if not hasattr(self, "_initial_log_dem"):
#         self._estimate_initial_dem()

#     results = {}
#     best_chi2 = np.inf
#     best_result = None
#     best_method = None

#     for method in methods:
#         print(f"\n>>> Trying method: {method}")
#         params = self._build_lmfit_parameters(n_knots=n_knots)
#         result = minimize(self._residuals, params, method=method, **kwargs)

#         # Compute DEM + chi square for this fit
#         # SELFNOTEJOY - output currently does not have units. unts=cm^5 * K^-1 Make this a test
#         dem = self._reconstruct_dem_from_knots(
#             result.params
#         )  # SELFNOTEJOY - here is the stamp to defining the DEM - triple check
#         dem_mid = 0.5 * (dem[:-1] + dem[1:])
#         R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
#         T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
#         I_fit = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

#         sigma = self.intensity_errors.to_value(u.DN / u.s)
#         residuals = (self._observed_intensities - I_fit) / sigma
#         chi2 = np.sum(residuals**2)

#         print(f"x square = {chi2:.3e}")

#         results[method] = (result, chi2)

#         if chi2 < best_chi2:
#             best_chi2 = chi2
#             best_result = result
#             best_method = method

#     print(f"\n>>> Best method: {best_method} with x square = {best_chi2:.3e}")

#     # Store outputs from the best fit
#     best_dem = self._reconstruct_dem_from_knots(best_result.params)
#     self.dem = best_dem
#     dem_mid = 0.5 * (best_dem[:-1] + best_dem[1:])
#     R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
#     T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
#     self.fitted_intensities = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)
#     sigma = self.intensity_errors.to_value(u.DN / u.s)
#     residuals = (self._observed_intensities - self.fitted_intensities) / sigma
#     self.chi2 = np.sum(residuals**2)
#     dof = len(self._observed_intensities) - len(best_result.params)
#     self.redchi2 = self.chi2 / max(dof, 1)

#     return best_result

# def run_monte_carlo(
#     self, n_runs=None, n_knots=6, method="least_squares", random_seed=None
# ):
#     from tqdm import tqdm  # add this at top of file

#     if n_runs is None:
#         n_runs = self._monte_carlo_runs
#     if n_runs <= 0:
#         raise ValueError("Monte Carlo runs disabled (n_runs=0).")

#     if random_seed is not None:
#         np.random.seed(random_seed)

#     sigma = self.intensity_errors.to_value(u.DN / u.s)
#     dem_ensemble = []

#     self._last_n_knots = n_knots

#     # --- progress bar
#     for i in tqdm(range(n_runs), desc="Monte Carlo DEM fits", unit="run"):
#         noisy_obs = self._observed_intensities + np.random.normal(0, sigma)
#         self._active_observed_intensities = noisy_obs
#         try:
#             params = self._build_lmfit_parameters(n_knots=n_knots)
#             result = minimize(self._residuals, params, method=method)
#         finally:
#             if hasattr(self, "_active_observed_intensities"):
#                 delattr(self, "_active_observed_intensities")

#         dem_i = self._reconstruct_dem_from_knots(result.params)
#         dem_ensemble.append(dem_i)

#     dem_ensemble = np.array(dem_ensemble)
#     self._dem_ensemble = dem_ensemble
#     self.dem_uncertainty = np.std(dem_ensemble, axis=0)
#     self.dem_median = np.median(dem_ensemble, axis=0)

#     return dem_ensemble

# # -------------------------------------------------------------------------------------------------------------------------------

# ################################################################################################################################
# ################################################################################################################################
# #############************************** Start of  error Summary **************************##########################
# ################## #####################

# def summary(self):
#     """
#     Print a comprehensive summary of the DEM solver setup,
#     including inputs, solver configuration, fit results,
#     Monte Carlo ensemble status, and available plotting helpers.
#     """
#     print("\nXRTpy DEM Iterative Setup Summary\n")
#     print("=" * 65)

#     # Filters & Observations
#     print(f" Filters:               {self.filter_names}")
#     print(f" Observed Intensities:  {self.observed_intensities}")
#     print(f" Number of channels:    {len(self._observed_intensities)}")

#     # Errors
#     print(f" Intensity Errors:      {self.intensity_errors}")
#     if self._intensity_errors is not None:
#         print(" Error model used:      User-provided")
#     else:
#         print(
#             f" Error model used:      Auto-estimated "
#             f"(obs * 0.03, min={self.min_observational_error.value} DN/s)"
#         )
#         print("   [IDL reference: xrt_dem_iterative2.pro]")

#     # Temperature grid
#     print(
#         f" Temperature grid:      logT {self.minimum_bound_temperature:.2f}–{self.maximum_bound_temperature:.2f}, step {self.logarithmic_temperature_step_size}"
#     )
#     print(f" Temp bins:             {len(self.logT)}")
#     print(f" dlogT:                 {self.dlogT:.3f}, dlnT: {self.dlnT:.3f}")

#     # Solver setup
#     print(f" Solver factor:         {self.normalization_factor:.1e}")
#     print(f" Monte Carlo runs:      {self.monte_carlo_runs or 'None'}")
#     print(f" Max Iterations:        {self.max_iterations}")
#     print(f" Knots (n_knots):       {getattr(self, '_last_n_knots', 'default=6')}")

#     if hasattr(self, "chi2"):
#         dof = len(self._observed_intensities) - len(
#             getattr(self, "_init_knot_params", [])
#         )
#         print(f"   χ²:                  {self.chi2:.4e} (dof={dof})")

#     # Responses
#     print(f" Response unit:         {self._response_unit}")
#     if hasattr(self, "_response_matrix"):
#         print(
#             f" Response matrix:       {self._response_matrix.shape} (filters × bins)"
#         )
#     else:
#         print(" Response matrix:       Not yet built")

#     # Fit results
#     if hasattr(self, "dem"):
#         print("\n Fit Results:")
#         print(f"   DEM bins:            {self.dem.shape}")
#         if hasattr(self, "chi2"):
#             print(f"   Chi²:                {self.chi2:.4e}")
#         if hasattr(self, "redchi2"):
#             print(f"   Reduced Chi²:        {self.redchi2:.4e}")
#         if hasattr(self, "fitted_intensities"):
#             print(f"   Fitted Intensities:  {self.fitted_intensities}")

#     # Monte Carlo results
#     if hasattr(self, "_dem_ensemble"):
#         print("\n Monte Carlo Ensemble:")
#         n_mc = len(self._dem_ensemble)
#         print(f"   Runs stored:         {n_mc}")
#         dem_stack = np.array(self._dem_ensemble)
#         med = np.median(dem_stack, axis=0)
#         spread = np.percentile(dem_stack, [16, 84], axis=0)
#         print("   DEM median (log10 cm^-5 K^-1):")
#         print(f"      First 5 bins:     {np.log10(med[:5]+1e-40)}")
#         print("   DEM 1σ spread (first bin):")
#         print(
#             f"      {np.log10(spread[0,0]+1e-40):.2f} – {np.log10(spread[1,0]+1e-40):.2f}"
#         )
#         print("   Reproducibility:     Run with random_seed for identical results")

#     if hasattr(self, "chi2"):
#         print(f"   Chi²:                {self.chi2:.4e}")
#     if hasattr(self, "redchi2"):
#         print(f"   Reduced Chi²:        {self.redchi2:.4e}")
#     if hasattr(self, "dof"):
#         print(f"   Degrees of Freedom:  {self.dof}")
#     if hasattr(self, "_iter_stats") and len(self._iter_stats["chisq"]) > 0:
#         print(f"   Iterations tracked:  {len(self._iter_stats['chisq'])}")
#         print(f"   Final Iter χ²:       {self._iter_stats['chisq'][-1]:.4e}")

#     # Plotting guidance
#     print("\n Plotting Options:")
#     if hasattr(self, "dem"):
#         print("   • plot_dem_results(results) → Quick plot from solve() dictionary")
#         print(
#             "   • plot_dem_uncertainty()   → Best-fit DEM + shaded ±1σ (if MC available)"
#         )
#         print(
#             "   • plot_idl_style()         → IDL-style view (best-fit + MC curves)"
#         )
#         print(
#             "   • plot_dem_with_median_bins() → Median + closest DEM (IDL style extension)"
#         )
#         print("   • plot_fit_residuals()     → Observed vs fitted intensities")
#         print("   • plot_iteration_stats()    ")

#     print("=" * 65)

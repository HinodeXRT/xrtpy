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
            self.response_temperatures, self.response_values
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

    def _estimate_initial_dem(self, cutoff: float = 1.0 / np.e) -> np.ndarray:
        """
        Estimate an initial DEM curve from observed intensities and responses.

        This follows the algorithm in IDL's `xrt_dem_iterative2.pro`, which uses
        response-peak inversion to generate a crude log10 DEM estimate per channel,
        then interpolates these estimates onto the solver's regular temperature grid.

        Parameters
        ----------
        cutoff : float, optional
            Fraction of the peak response to use for defining the "good" window
            around each filter's peak. Default is 1/e ≈ 0.3679.

        Returns
        -------
        est_log_dem_on_grid : ndarray
            Array of shape (n_temperatures,) giving the initial DEM estimate
            on `self.logT`. Values are log10(DEM) in [cm^-5 K^-1].
            This can be used to seed the solver.

        Notes
        -----
        - Units:
            * Observed intensities: [DN s^-1 pix^-1]
            * Response: [DN s^-1 pix^-1 cm^5]
            * DEM(T): [cm^-5 K^-1]
        - For each filter:
            1. Locate the peak of its response.
            2. Define a window where response > cutoff * peak.
            3. Compute the denominator integral: sum( T * R * dlnT ).
            4. Estimate DEM_peak = I_obs / denom.
            5. Store log10(DEM_peak) at the peak logT.
        - Duplicate peak logTs are merged by averaging.
        - If fewer than 2 valid points are found, falls back to a flat guess
        (log10 DEM = 22 everywhere).
        """
        if not hasattr(self, "logT"):
            raise AttributeError(
                "Temperature grid missing. Call create_logT_grid() first."
            )
        if not hasattr(self, "_response_matrix"):
            raise AttributeError(
                "Response matrix missing. Call _interpolate_responses_to_grid() first."
            )

        # Storage for peak locations and DEM estimates
        t_peaks = []
        log_dem_estimates = []

        # Loop over each filter
        for i, (T_orig, R_orig, I_obs) in enumerate(
            zip(
                self.response_temperatures,
                self.response_values,
                self._observed_intensities,
            )
        ):
            logT_orig = np.log10(T_orig.to_value(u.K))
            R_vals = R_orig.to_value((u.DN / u.s / u.pix) * u.cm**5)

            if I_obs <= 0 or np.all(R_vals <= 0):
                continue  # skip unusable channel

            # 1. Peak location
            max_idx = np.argmax(R_vals)
            peak_val = R_vals[max_idx]
            t_peak = (
                np.round(logT_orig[max_idx] / self._logarithmic_temperature_step_size)
                * self._logarithmic_temperature_step_size
            )

            # 2. Good window (where R > cutoff * peak)
            good = np.where(R_vals > peak_val * cutoff)[0]
            if len(good) < 2:
                continue

            # 3. Compute denominator integral: sum(T * R * dlnT)
            T_good = 10.0 ** logT_orig[good]
            R_good = R_vals[good]
            dlogT_native = np.diff(logT_orig).mean()
            dlnT_native = np.log(10.0) * dlogT_native
            denom = np.sum(T_good * R_good * dlnT_native)

            if denom <= 0:
                continue

            # 4. DEM estimate at peak
            dem_peak = I_obs / denom  # [cm^-5 K^-1]
            if dem_peak <= 0 or not np.isfinite(dem_peak):
                continue

            log_dem_est = np.log10(dem_peak)
            t_peaks.append(t_peak)
            log_dem_estimates.append(log_dem_est)

        # 5. Handle duplicates: average log10 DEM at same t_peak
        if len(t_peaks) == 0:
            # Fallback: flat guess (IDL style)
            est_log_dem_on_grid = np.ones_like(self.logT) * 22.0
            self._initial_log_dem = est_log_dem_on_grid
            return est_log_dem_on_grid

        uniq_t = {}
        for t, dem_val in zip(t_peaks, log_dem_estimates):
            if t in uniq_t:
                uniq_t[t].append(dem_val)
            else:
                uniq_t[t] = [dem_val]
        t_peaks_uniq = np.array(sorted(uniq_t.keys()))
        log_dem_uniq = np.array([np.mean(uniq_t[t]) for t in t_peaks_uniq])

        if len(t_peaks_uniq) < 2:
            # Not enough points > flat guess
            est_log_dem_on_grid = np.ones_like(self.logT) * 22.0
            self._initial_log_dem = est_log_dem_on_grid
            return est_log_dem_on_grid

        # 6. Interpolate sparse estimates onto the solver's grid
        interp_func = interp1d(
            t_peaks_uniq,
            log_dem_uniq,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        est_log_dem_on_grid = interp_func(self.logT)

        # Store for later use
        self._initial_log_dem = est_log_dem_on_grid

        return est_log_dem_on_grid

    def _build_lmfit_parameters(self, n_knots: int = 6):
        """
        Build lmfit.Parameters for the DEM spline knots.

        Parameters
        ----------
        n_knots : int, optional
            Number of spline knots across the logT grid. Default = 6.

        Returns
        -------
        params : lmfit.Parameters
            Parameters object containing log10(DEM/normalization_factor) values at knot points.
            Each parameter is named "knot_i" where i = 0..n_knots-1.

        Notes
        -----
        - IDL's `xrt_dem_iterative2.pro` seeds its fit by taking DEM estimates
        at peak response temperatures and spreading them across the grid.
        - Here, we select evenly spaced knots across the solver's logT range.
        - The stored value at each knot is:
            log10(DEM / normalization_factor)
        where `normalization_factor` is typically 1e17.
        - Bounds can be applied to prevent extreme DEM excursions if desired.
        """
        if not hasattr(self, "_initial_log_dem"):
            raise AttributeError(
                "Initial DEM not available. Run _estimate_initial_dem() first."
            )

        from lmfit import Parameters

        # Choose evenly spaced knot positions across logT range
        knot_positions = np.linspace(
            self._minimum_bound_temperature, self._maximum_bound_temperature, n_knots
        )
        self._knot_positions = knot_positions  # store for later reconstruction

        # Interpolate initial DEM estimate at these knot positions
        interp_func = interp1d(
            self.logT,
            self._initial_log_dem,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        init_log_dem_at_knots = interp_func(knot_positions)

        # Convert to log10(DEM/normalization_factor)
        init_scaled = init_log_dem_at_knots - np.log10(self._normalization_factor)

        # Build lmfit Parameters
        params = Parameters()
        for i, val in enumerate(init_scaled):
            params.add(
                name=f"knot_{i}",
                value=val,
                min=-10,  # optional bounds: avoid absurdly low
                max=50,  # optional bounds: avoid absurdly high
                vary=True,
            )

        self._init_knot_params = params
        return params

    def _reconstruct_dem_from_knots(self, params) -> np.ndarray:
        """
        Reconstruct the DEM curve on the solver's logT grid from spline knot parameters.

        Parameters
        ----------
        params : lmfit.Parameters
            Knot parameters where each value is log10(DEM / normalization_factor).

        Returns
        -------
        dem_grid : ndarray
            DEM values on `self.logT` grid in linear space [cm^-5 K^-1].

        Notes
        -----
        - Knot positions are stored in `self._knot_positions` when
        `_build_lmfit_parameters()` is called.
        - The stored parameter values are log10(DEM/normalization_factor).
        - Conversion back to DEM:
            DEM = normalization_factor * 10^(interp(log10 DEM/normalization_factor))
        - Interpolation is linear in log space, as in IDL's `xrt_dem_iterative2.pro`.
        """
        if not hasattr(self, "_knot_positions"):
            raise AttributeError(
                "Knot positions not found. Run _build_lmfit_parameters() first."
            )

        # Extract knot values from parameters (log10(DEM/normalization_factor))
        knot_vals = np.array(
            [params[f"knot_{i}"].value for i in range(len(self._knot_positions))]
        )

        # Interpolate across solver grid in log space
        interp_func = interp1d(
            self._knot_positions,
            knot_vals,
            kind="linear",
            bounds_error=False,
            fill_value="extrapolate",
        )
        log_dem_scaled = interp_func(self.logT)

        # Convert back to DEM [cm^-5 K^-1]
        dem_grid = self._normalization_factor * (
            10.0**log_dem_scaled
        )  ## dem_grid now back in physical units

        return dem_grid

    # self._iteration_chi2 = []

    # def _residuals(self, params) -> np.ndarray:
    #     """
    #     Residuals function for DEM fitting.

    #     Parameters
    #     ----------
    #     params : lmfit.Parameters
    #         Knot parameters, each storing log10(DEM / normalization_factor).

    #     Returns
    #     -------
    #     residuals : ndarray
    #         Vector of normalized residuals for each observed channel:
    #         (I_obs - I_calc) / sigma
    #         Shape = (n_filters,)

    #     Notes
    #     -----
    #     - This is the core of the DEM solver. It reconstructs the DEM curve
    #     from spline knot parameters, computes modeled intensities by
    #     integrating DEM * Response over temperature, and returns the
    #     residuals relative to observations.
    #     - Integration is done using midpoint trapezoid approximation:
    #         I_calc[i] = sum_j DEM_mid[j] * R_mid[i,j] * T_mid[j] * dlnT
    #     """
    #     # 1. Reconstruct DEM on the grid
    #     dem_grid = self._reconstruct_dem_from_knots(params)  # [cm^-5 K^-1]

    #     # 2. Prepare midpoint arrays
    #     dem_mid = 0.5 * (dem_grid[:-1] + dem_grid[1:])
    #     R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
    #     T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)

    #     # 3. Compute modeled intensities
    #     # Shape: (n_filters,)
    #     I_calc = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

    #     # 4. Residuals: normalize by observational errors
    #     # sigma = self.intensity_errors.to_value(u.DN / u.s)  # ensure numeric
    #     # residuals = (self._observed_intensities - I_calc) / sigma
    #     # Use MC-perturbed intensities if present; otherwise the originals

    #     y_obs = getattr(self, "_active_observed_intensities", self._observed_intensities)
    #     sigma = self.intensity_errors.to_value(u.DN / u.s)  # numeric
    #     residuals = (y_obs - I_calc) / sigma

    #     residuals = (y_obs - I_calc) / sigma
    #     chi2_val = np.sum(residuals**2)

    #     # Log χ² per iteration
    #     if not hasattr(self, "_iteration_chi2"):
    #         self._iteration_chi2 = []
    #     self._iteration_chi2.append(chi2_val)

    #     return residuals
    def _residuals(self, params) -> np.ndarray:
        """
        Residuals function for DEM fitting.

        Returns
        -------
        residuals : ndarray
            (I_obs - I_calc) / sigma, one per filter.
        """
        # 1. Reconstruct DEM on the grid
        dem_grid = self._reconstruct_dem_from_knots(params)  # [cm^-5 K^-1]

        # 2. Midpoint integration setup
        dem_mid = 0.5 * (dem_grid[:-1] + dem_grid[1:])
        R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
        T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)

        # 3. Modeled intensities
        I_calc = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

        # # 4. Residuals: normalize by observational errors
        # y_obs = getattr(
        #     self, "_active_observed_intensities", self._observed_intensities
        # )
        # sigma = self.intensity_errors.to_value(u.DN / u.s)
        # residuals = (y_obs - I_calc) / sigma

        # # 5. Track χ² per iteration
        # chi2_val = np.sum(residuals**2)
        # if not hasattr(self, "_iteration_chi2"):
        #     self._iteration_chi2 = []
        # self._iteration_chi2.append(chi2_val)

        # return residuals

        # 4. Use either base or MC-perturbed observed intensities (physical units)
        y_obs_phys = getattr(
            self, "_active_observed_intensities", self._observed_intensities
        )
        sigma_phys = self.intensity_errors.to_value(u.DN / u.s)

        # 5. Apply normalization_factor in an IDL-like way:
        #    scale data, model, and errors by the same factor.
        #    (This keeps residuals numerically identical, but keeps
        #     internal numbers closer to order unity if desired.)
        nf = self._normalization_factor

        y_scaled = y_obs_phys / nf
        I_calc_scaled = I_calc / nf
        sigma_scaled = sigma_phys / nf

        residuals = (y_scaled - I_calc_scaled) / sigma_scaled

        # 6. Track χ² per iteration (in scaled space — same χ² as unscaled)
        chi2_val = np.sum(residuals**2)
        if not hasattr(self, "_iteration_chi2"):
            self._iteration_chi2 = []
        self._iteration_chi2.append(chi2_val)

        return residuals

    # def fit_dem(self, n_knots: int = 6, method: str = "least_squares", **kwargs):
    #     """
    #     Fit the DEM using lmfit to minimize residuals.

    #     Parameters
    #     ----------
    #     n_knots : int, optional
    #         Number of spline knots across the logT grid. Default = 6.
    #     method : str, optional
    #         Minimization method passed to `lmfit.minimize`.
    #         Common choices: "least_squares", "leastsq", "nelder".
    #         Default = "least_squares".
    #     **kwargs : dict
    #         Additional keyword arguments forwarded to `lmfit.minimize`.

    #     Returns
    #     -------
    #     result : lmfit.MinimizerResult
    #         The lmfit result object containing fit information.

    #     Side Effects
    #     ------------
    #     On successful fit, stores:
    #     - self.dem : ndarray
    #         Best-fit DEM(T) on self.logT [cm^-5 K^-1].
    #     - self.fitted_intensities : ndarray
    #         Modeled intensities [DN/s/pix] for each filter.
    #     - self.chi2 : float
    #         Chi-squared (sum of squared residuals).
    #     - self.redchi2 : float
    #         Reduced chi-squared, normalized by (Nobs - Nparams).

    #     Notes
    #     -----
    #     This method automatically builds the logT grid, interpolates
    #     responses, and estimates an initial DEM if not already done.
    #     """
    #     from lmfit import minimize

    #     # --- Auto-prepare prerequisites ---
    #     if not hasattr(self, "logT") or not hasattr(self, "T"):
    #         self.create_logT_grid()

    #     if not hasattr(self, "_response_matrix"):
    #         self._interpolate_responses_to_grid()

    #     if not hasattr(self, "_initial_log_dem"):
    #         self._estimate_initial_dem()

    #     self._last_n_knots = n_knots #Used for print in the summary function

    #     # 1. Build initial knot parameters
    #     params = self._build_lmfit_parameters(n_knots=n_knots)

    #     # 2. Run minimization
    #     result = minimize(self._residuals, params, method=method, **kwargs)

    #     # 3. On success, reconstruct DEM and fitted intensities
    #     best_dem = self._reconstruct_dem_from_knots(result.params)
    #     self.dem = best_dem  # [cm^-5 K^-1]

    #     # Compute fitted intensities using midpoint integration
    #     dem_mid = 0.5 * (best_dem[:-1] + best_dem[1:])
    #     R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
    #     T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
    #     I_fit = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

    #     self.fitted_intensities = I_fit  # [DN/s/pix]
    #     sigma = self.intensity_errors.to_value(u.DN / u.s)
    #     residuals = (self._observed_intensities - I_fit) / sigma

    #     # Chi-squared metrics
    #     self.chi2 = np.sum(residuals**2)
    #     dof = len(self._observed_intensities) - len(result.params)
    #     self.redchi2 = self.chi2 / max(dof, 1)

    #     return result

    def fit_dem(self, n_knots: int = 6, method: str = "least_squares", **kwargs):
        """
        Fit the DEM using lmfit to minimize residuals.
        Tracks chi² per iteration (like IDL's XRT_ITER_DEMSTAT).
        """
        from lmfit import minimize

        # --- Auto-prepare prerequisites ---
        if not hasattr(self, "logT") or not hasattr(self, "T"):
            self.create_logT_grid()
        if not hasattr(self, "_response_matrix"):
            self._interpolate_responses_to_grid()
        if not hasattr(self, "_initial_log_dem"):
            self._estimate_initial_dem()

        self._last_n_knots = n_knots  # for summary()

        # Storage for iteration statistics
        self._iter_stats = {"chisq": [], "iteration": []}

        def _callback(params, iter, resid, *args, **kwargs):
            # Compute chi² at this iteration
            chi2 = np.sum(resid**2)
            self._iter_stats["chisq"].append(chi2)
            self._iter_stats["iteration"].append(iter)

        # 1. Build initial knot parameters
        params = self._build_lmfit_parameters(n_knots=n_knots)

        # 2. Run minimization
        result = minimize(
            self._residuals,
            params,
            method=method,
            iter_cb=_callback,  # <-- track stats
            max_nfev=self.max_iterations,
            **kwargs,
        )

        # 3. On success, reconstruct DEM + fitted intensities
        best_dem = self._reconstruct_dem_from_knots(result.params)
        self.dem = best_dem

        dem_mid = 0.5 * (best_dem[:-1] + best_dem[1:])
        R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
        T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
        I_fit = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

        self.fitted_intensities = I_fit
        sigma = self.intensity_errors.to_value(u.DN / u.s)
        residuals = (self._observed_intensities - I_fit) / sigma

        # Chi^2 metrics
        self.chi2 = np.sum(residuals**2)
        dof = len(self._observed_intensities) - len(result.params)
        self.redchi2 = self.chi2 / max(dof, 1)
        self.dof = dof  # save for summary

        return result

    def fit_with_multiple_methods(
        self, methods=("leastsq", "least_squares", "nelder"), n_knots: int = 6, **kwargs
    ):
        """
        Try multiple lmfit minimization methods and pick the best χ².

        Parameters
        ----------
        methods : tuple of str, optional
            Minimization methods to test. Default = ("leastsq", "least_squares", "nelder").
        n_knots : int, optional
            Number of spline knots for DEM fit. Default = 6.
        **kwargs : dict
            Extra arguments passed to `lmfit.minimize`.

        Returns
        -------
        best_result : lmfit.MinimizerResult
            Result from the method with lowest chi².
        """
        from lmfit import minimize

        if not hasattr(self, "_initial_log_dem"):
            self._estimate_initial_dem()

        results = {}
        best_chi2 = np.inf
        best_result = None
        best_method = None

        for method in methods:
            print(f"\n>>> Trying method: {method}")
            params = self._build_lmfit_parameters(n_knots=n_knots)
            result = minimize(self._residuals, params, method=method, **kwargs)

            # Compute DEM + chi square for this fit
            # SELFNOTEJOY - output currently does not have units. unts=cm^5 * K^-1 Make this a test
            dem = self._reconstruct_dem_from_knots(
                result.params
            )  # SELFNOTEJOY - here is the stamp to defining the DEM - triple check
            dem_mid = 0.5 * (dem[:-1] + dem[1:])
            R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
            T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
            I_fit = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)

            sigma = self.intensity_errors.to_value(u.DN / u.s)
            residuals = (self._observed_intensities - I_fit) / sigma
            chi2 = np.sum(residuals**2)

            print(f"x square = {chi2:.3e}")

            results[method] = (result, chi2)

            if chi2 < best_chi2:
                best_chi2 = chi2
                best_result = result
                best_method = method

        print(f"\n>>> Best method: {best_method} with x square = {best_chi2:.3e}")

        # Store outputs from the best fit
        best_dem = self._reconstruct_dem_from_knots(best_result.params)
        self.dem = best_dem
        dem_mid = 0.5 * (best_dem[:-1] + best_dem[1:])
        R_mid = 0.5 * (self._response_matrix[:, :-1] + self._response_matrix[:, 1:])
        T_mid = 0.5 * (self.T[:-1] + self.T[1:]).to_value(u.K)
        self.fitted_intensities = np.sum(R_mid * dem_mid * T_mid * self.dlnT, axis=1)
        sigma = self.intensity_errors.to_value(u.DN / u.s)
        residuals = (self._observed_intensities - self.fitted_intensities) / sigma
        self.chi2 = np.sum(residuals**2)
        dof = len(self._observed_intensities) - len(best_result.params)
        self.redchi2 = self.chi2 / max(dof, 1)

        return best_result

    # def run_monte_carlo(self, n_runs=None, n_knots=6, method="least_squares", random_seed=None):
    #     if random_seed is not None:
    #         np.random.seed(random_seed)
    #     """
    #     Run Monte Carlo DEM fits to estimate uncertainties and store full ensemble.

    #     Returns
    #     -------
    #     dem_ensemble : ndarray
    #         Shape (n_runs, n_temperatures) array of DEM solutions.
    #     """
    #     from lmfit import minimize

    #     if n_runs is None:
    #         n_runs = self._monte_carlo_runs
    #     if n_runs <= 0:
    #         raise ValueError("Monte Carlo runs disabled (n_runs=0).")

    #     sigma = self.intensity_errors.to_value(u.DN / u.s)
    #     dem_ensemble = []

    #     self._last_n_knots = n_knots #Used for print in the summary function

    #     for i in range(n_runs):
    #         noisy_obs = self._observed_intensities + np.random.normal(0, sigma)
    #         #print(f"Given intensities: {noisy_obs}")
    #         self._observed_intensities_mc = noisy_obs  # temp override

    #         # params = self._build_lmfit_parameters(n_knots=n_knots)  #Older Version  Sept 18
    #         # result = minimize(lambda p: self._residuals(p), params, method=method) #Older Version Sept 18
    #         params = self._build_lmfit_parameters(n_knots=n_knots)
    #         # Activate noisy intensities for this run
    #         self._active_observed_intensities = noisy_obs
    #         try:
    #             result = minimize(self._residuals, params, method=method)
    #         finally:
    #             # Always restore (so the main dataset isn’t polluted)
    #             if hasattr(self, "_active_observed_intensities"):
    #                 delattr(self, "_active_observed_intensities")

    #         dem_i = self._reconstruct_dem_from_knots(result.params)
    #         dem_ensemble.append(dem_i)

    #     dem_ensemble = np.array(dem_ensemble)

    #     # Store ensemble + uncertainty
    #     self._dem_ensemble = dem_ensemble
    #     self.dem_uncertainty = np.std(dem_ensemble, axis=0)
    #     self.dem_median = np.median(dem_ensemble, axis=0)

    #     return dem_ensemble

    def run_monte_carlo(
        self, n_runs=None, n_knots=6, method="least_squares", random_seed=None
    ):
        import numpy as np
        from lmfit import minimize
        from tqdm import tqdm  # add this at top of file

        if n_runs is None:
            n_runs = self._monte_carlo_runs
        if n_runs <= 0:
            raise ValueError("Monte Carlo runs disabled (n_runs=0).")

        if random_seed is not None:
            np.random.seed(random_seed)

        sigma = self.intensity_errors.to_value(u.DN / u.s)
        dem_ensemble = []

        self._last_n_knots = n_knots

        # --- progress bar
        for i in tqdm(range(n_runs), desc="Monte Carlo DEM fits", unit="run"):
            noisy_obs = self._observed_intensities + np.random.normal(0, sigma)
            self._active_observed_intensities = noisy_obs
            try:
                params = self._build_lmfit_parameters(n_knots=n_knots)
                result = minimize(self._residuals, params, method=method)
            finally:
                if hasattr(self, "_active_observed_intensities"):
                    delattr(self, "_active_observed_intensities")

            dem_i = self._reconstruct_dem_from_knots(result.params)
            dem_ensemble.append(dem_i)

        dem_ensemble = np.array(dem_ensemble)
        self._dem_ensemble = dem_ensemble
        self.dem_uncertainty = np.std(dem_ensemble, axis=0)
        self.dem_median = np.median(dem_ensemble, axis=0)

        return dem_ensemble

    def solve(
        self, n_knots: int = 6, method: str = "least_squares", run_mc: bool = True
    ):
        """
        Run the full DEM solver, IDL-style.

        This orchestrates:
        1. Build temperature grid.
        2. Interpolate responses (response matrix).
        3. Estimate initial DEM.
        4. Fit DEM with lmfit.
        5. Optionally run Monte Carlo ensemble.

        Parameters
        ----------
        n_knots : int, optional
            Number of spline knots across logT. Default = 6.
        method : str, optional
            Minimization method for `lmfit.minimize`. Default = "least_squares".
        run_mc : bool, optional
            Whether to run Monte Carlo simulations (using self.monte_carlo_runs).
            Default = True.

        Returns
        -------
        results : dict
            Dictionary of solver outputs:
            - "temperature" : log10(T) grid
            - "dem"         : best-fit DEM [cm^-5 K^-1]
            - "dem_err"     : DEM uncertainty (if MC runs > 0)
            - "ifit"        : fitted intensities [DN/s/pix]
            - "chi2"        : χ²
            - "redchi2"     : reduced χ²
        """
        # IDL-STYLE NOSOLVE CHECk
        # IDL behavior: if all observed intensities are zero (or non-positive),
        # the DEM is trivially zero. Skip solving and return immediately.
        if np.all(self._observed_intensities <= 0):  # == 0
            warnings.warn(
                "\n\n All observed intensities are zero or non-positive. "
                "DEM cannot be solved. Returning zero DEM and zero fitted intensities. \n\n"
            )

            # Ensure grid exists (IDL also returns logT_out even for nosolve)
            if not hasattr(self, "logT"):
                self.create_logT_grid()

            self.dem = np.zeros_like(self.logT)
            self.fitted_intensities = np.zeros_like(self._observed_intensities)
            self.chi2 = 0.0
            self.redchi2 = 0.0
            return self

        # 1. Ensure grid & responses
        self.create_logT_grid()
        self._interpolate_responses_to_grid()

        # 2. Estimate initial DEM
        self._estimate_initial_dem()

        # 3. Fit DEM
        result = self.fit_dem(n_knots=n_knots, method=method)

        # 4. Monte Carlo (optional)
        if run_mc and self.monte_carlo_runs > 0:
            self.run_monte_carlo(
                n_runs=self.monte_carlo_runs, n_knots=n_knots, method=method
            )

        # # 5. Bundle results
        # return {
        #     "temperature": self.logT,
        #     "dem": self.dem,
        #     "dem_err": getattr(self, "dem_uncertainty", None),
        #     "ifit": self.fitted_intensities,
        #     "chi2": getattr(self, "chi2", None),
        #     "redchi2": getattr(self, "redchi2", None),
        #     "solver": self,
        # }
        return self

    def to_dict(self):
        """Return solver outputs as a dictionary."""
        return {
            "temperature": self.logT,
            "dem": getattr(self, "dem", None),
            "dem_err": getattr(self, "dem_uncertainty", None),
            "ifit": getattr(self, "fitted_intensities", None),
            "chi2": getattr(self, "chi2", None),
            "redchi2": getattr(self, "redchi2", None),
        }

    def summary(self):
        """
        Print a comprehensive summary of the DEM solver setup,
        including inputs, solver configuration, fit results,
        Monte Carlo ensemble status, and available plotting helpers.
        """
        print("\nXRTpy DEM Iterative Setup Summary\n")
        print("=" * 65)

        # Filters & Observations
        print(f" Filters:               {self.filter_names}")
        print(f" Observed Intensities:  {self.observed_intensities}")
        print(f" Number of channels:    {len(self._observed_intensities)}")

        # Errors
        print(f" Intensity Errors:      {self.intensity_errors}")
        if self._intensity_errors is not None:
            print(" Error model used:      User-provided")
        else:
            print(
                f" Error model used:      Auto-estimated "
                f"(obs * 0.03, min={self.min_observational_error.value} DN/s)"
            )
            print("   [IDL reference: xrt_dem_iterative2.pro]")

        # Temperature grid
        print(
            f" Temperature grid:      logT {self.minimum_bound_temperature:.2f}–{self.maximum_bound_temperature:.2f}, step {self.logarithmic_temperature_step_size}"
        )
        print(f" Temp bins:             {len(self.logT)}")
        print(f" dlogT:                 {self.dlogT:.3f}, dlnT: {self.dlnT:.3f}")

        # Solver setup
        print(f" Solver factor:         {self.normalization_factor:.1e}")
        print(f" Monte Carlo runs:      {self.monte_carlo_runs or 'None'}")
        print(f" Max Iterations:        {self.max_iterations}")
        print(f" Knots (n_knots):       {getattr(self, '_last_n_knots', 'default=6')}")

        if hasattr(self, "chi2"):
            dof = len(self._observed_intensities) - len(
                getattr(self, "_init_knot_params", [])
            )
            print(f"   χ²:                  {self.chi2:.4e} (dof={dof})")

        # Responses
        print(f" Response unit:         {self._response_unit}")
        if hasattr(self, "_response_matrix"):
            print(
                f" Response matrix:       {self._response_matrix.shape} (filters × bins)"
            )
        else:
            print(" Response matrix:       Not yet built")

        # Fit results
        if hasattr(self, "dem"):
            print("\n Fit Results:")
            print(f"   DEM bins:            {self.dem.shape}")
            if hasattr(self, "chi2"):
                print(f"   Chi²:                {self.chi2:.4e}")
            if hasattr(self, "redchi2"):
                print(f"   Reduced Chi²:        {self.redchi2:.4e}")
            if hasattr(self, "fitted_intensities"):
                print(f"   Fitted Intensities:  {self.fitted_intensities}")

        # Monte Carlo results
        if hasattr(self, "_dem_ensemble"):
            print("\n Monte Carlo Ensemble:")
            n_mc = len(self._dem_ensemble)
            print(f"   Runs stored:         {n_mc}")
            dem_stack = np.array(self._dem_ensemble)
            med = np.median(dem_stack, axis=0)
            spread = np.percentile(dem_stack, [16, 84], axis=0)
            print("   DEM median (log10 cm^-5 K^-1):")
            print(f"      First 5 bins:     {np.log10(med[:5]+1e-40)}")
            print("   DEM 1σ spread (first bin):")
            print(
                f"      {np.log10(spread[0,0]+1e-40):.2f} – {np.log10(spread[1,0]+1e-40):.2f}"
            )
            print("   Reproducibility:     Run with random_seed for identical results")

        if hasattr(self, "chi2"):
            print(f"   Chi²:                {self.chi2:.4e}")
        if hasattr(self, "redchi2"):
            print(f"   Reduced Chi²:        {self.redchi2:.4e}")
        if hasattr(self, "dof"):
            print(f"   Degrees of Freedom:  {self.dof}")
        if hasattr(self, "_iter_stats") and len(self._iter_stats["chisq"]) > 0:
            print(f"   Iterations tracked:  {len(self._iter_stats['chisq'])}")
            print(f"   Final Iter χ²:       {self._iter_stats['chisq'][-1]:.4e}")

        # Plotting guidance
        print("\n Plotting Options:")
        if hasattr(self, "dem"):
            print("   • plot_dem_results(results) → Quick plot from solve() dictionary")
            print(
                "   • plot_dem_uncertainty()   → Best-fit DEM + shaded ±1σ (if MC available)"
            )
            print(
                "   • plot_idl_style()         → IDL-style view (best-fit + MC curves)"
            )
            print(
                "   • plot_dem_with_median_bins() → Median + closest DEM (IDL style extension)"
            )
            print("   • plot_fit_residuals()     → Observed vs fitted intensities")
            print("   • plot_iteration_stats()    ")

        print("=" * 65)


# Attach plotting functions from plotting.py to the class
XRTDEMIterative.plot_dem_results = dem_plotting.plot_dem_results
XRTDEMIterative.plot_dem_uncertainty = dem_plotting.plot_dem_uncertainty
XRTDEMIterative.plot_idl_style = dem_plotting.plot_idl_style
XRTDEMIterative.plot_fit_residuals = dem_plotting.plot_fit_residuals
XRTDEMIterative.plot_dem_with_median_bins = dem_plotting.plot_dem_with_median_bins
XRTDEMIterative.plot_iteration_stats = dem_plotting.plot_iteration_stats


# NOTEFROMJOYTOJOY
# Missin outputs
# 1
# I'm missing BASE_OBS- It is a 2D array of intensity values for all DEM runs (base + MC)
# the set of observed intensities used in each DEM solution, including:
# Column 0 > the original observed intensities (your real data)
# Columns 1..MC_ITER > each Monte-Carlo–perturbed intensity vector

# 2
# MOD_OBS - Model intensities predicted by the DEM for each run.

# 3
# CHISQ- [1+MC_iter]
# chi-squre for each DEM solution
# Computed via mpdemfunct residuals squared and summed

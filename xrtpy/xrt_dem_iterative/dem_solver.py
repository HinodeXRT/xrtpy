__all__ = [
    "XRTDEMIterative",
]

import warnings

import astropy.units as u
import numpy as np
from lmfit import Parameters, minimize
from scipy.interpolate import interp1d, CubicSpline

from xrtpy.util.filters import validate_and_format_filters


class XRTDEMIterative:
    """
    Estimate the differential emission measure (DEM) from Hinode/XRT data
    using the iterative spline-based method.

    Parameters
    ----------
    observed_channel : str or list of str
        Filter names used in the observation (e.g., 'Al-mesh', 'Be-thin').
    observed_intensities : array-like
        Observed intensities in DN/s/pix for each channel.
    temperature_responses : list
        List of `TemperatureResponseFundamental` objects matching the filters.
    intensity_errors : array-like, optional
        Intensity uncertainties. If None, will use a model-based estimate.
    min_T : float
        Minimum log10 temperature (default: 5.5).
    max_T : float
        Maximum log10 temperature (default: 8.0).
    dT : float
        Step size in log10 temperature space (default: 0.1).
    min_error : float
        Minimum absolute intensity error (default: 2.0 DN/s/pix).
    relative_error : float
        Relative error for model-based uncertainty estimate (default: 0.03).
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
        min_T=5.5,
        max_T=8.0,
        dT=0.1,
        min_error=2.0,
        relative_error=0.03,
        monte_carlo_runs=0,
        max_iterations=2000,
        solv_factor=1e21,
    ):
        """
        Args:
            observed_channel (_type_): _description_
            observed_intensities (_type_): _description_
            temperature_responses (_type_): _description_
            intensity_errors (_type_, optional): _description_. Defaults to None.
            min_T (float, optional): _description_. Defaults to 5.5.
            max_T (float, optional): _description_. Defaults to 8.0.
            dT (float, optional): _description_. Defaults to 0.1.
            min_error (float, optional): _description_. Defaults to 2.0.
            relative_error (float, optional): _description_. Defaults to 0.03.

        Notes
        -----
        - All input lists (`observed_channel`, `observed_intensities`, and `temperature_responses`)
        must be the same length. Each entry should correspond to one filter.

        - The temperature grid range (`min_T`, `max_T`) must lie entirely within the
        response temperature ranges for **all** filters provided.

        - If `intensity_errors` is not provided, a model-based error estimate is used:
        max(relative_error * observed_intensity, min_error), as in the IDL original.

        - Default XRT filter names include:
        {'Al-mesh', 'Al-poly', 'C-poly', 'Ti-poly', 'Be-thin', 'Be-med', 'Al-med', 'Al-thick', 'Be-thick',
        'Al-poly/Al-mesh', 'Al-poly/Ti-poly', 'Al-poly/Al-thick', 'Al-poly/Be-thick'}
        """
        
        # Validate and store filter names
        if observed_channel is None or (hasattr(observed_channel, "__len__") and len(observed_channel) == 0):
            raise ValueError("`observed_channel` is required and cannot be empty.")
        self.observed_channel = validate_and_format_filters(observed_channel)


        # Store intensity and error arrays
        self._observed_intensities = np.asarray(observed_intensities, dtype=float)

        if observed_intensities is None or len(observed_intensities) == 0:
            raise ValueError("`observed_intensities` is required and cannot be empty.")
        
        if not np.all(np.isfinite(self._observed_intensities)):
            raise ValueError("`observed_intensities` must be finite numbers.")

        # Errors
        if intensity_errors is not None:
            self._intensity_errors = np.asarray(intensity_errors, dtype=float)
            if self._intensity_errors.shape != self._observed_intensities.shape:
                raise ValueError("Length of intensity_errors must match observed_intensities.")
            if not np.all(np.isfinite(self._intensity_errors)) or np.any(self._intensity_errors < 0):
                raise ValueError("`intensity_errors` must be finite and >= 0.")
        else:
            self._intensity_errors = None # Will be computed later




        # Store temperature grid parameters
        self._min_T = float(min_T)
        self._max_T = float(max_T)
        if not (self._min_T < self._max_T):
            raise ValueError("min_T must be < max_T.")
        if n_pts < 4:
            raise ValueError("Temperature grid must have at least 4 points.")

        self._dT = float(dT)

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

        # Check dT is positive
        if self._dT <= 0:
            raise ValueError("dT must be a positive scalar.")

        # Store temperature response objects
        self.responses = temperature_responses

        if temperature_responses is None or len(temperature_responses) == 0:
            raise ValueError("`temperature_responses` is required and cannot be empty.")

        # Validate that the temperature grid falls within the responses
        for r in self.responses:
            logT_grid = np.log10(r.temperature.value)
            if not (self._min_T >= logT_grid.min() and self._max_T <= logT_grid.max()):
                raise ValueError(
                    f"The specified temperature range [{min_T}, {max_T}] is outside the bounds of one or more filter response grids.\n"
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

        self.logT = np.arange(self._min_T, self._max_T + self._dT / 2, self._dT)

        # Store error model parameters
        self._min_error = float(min_error)
        self._relative_error = float(relative_error)

        # Validate and store intensity errors
        if intensity_errors is not None:
            self._intensity_errors = np.asarray(intensity_errors, dtype=float)
            if self._intensity_errors.shape != self._observed_intensities.shape:
                raise ValueError(
                    "Length of intensity_errors must match observed_intensities."
                )
        else:
            self._intensity_errors = None

        try:
            self._solv_factor = float(solv_factor)
            if self._solv_factor <= 0:
                raise ValueError("solv_factor must be a positive number.")
        except Exception as e:
            raise ValueError(f"Invalid solv_factor: {e}")

    def __repr__(self):
        return f"<XRTDEMIterative(filters={self.filter_names}, logT={self.min_T}–{self.max_T}, dT={self.dT})>"

    # @property  #Removed if not used
    # def name(self) -> str:
    #     """
    #     The XRT filter channel name, standardized (e.g. "Al-mesh").
    #     """
    #     return self._name

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
            Intensities in DN/s for each filter channel.
        """
        return self._observed_intensities * (u.DN / u.s)

    @property
    def filter_names(self):
        """
        Returns a list of filter names from the temperature responses.
        """
        return [r.filter_name for r in self.responses]

    @property
    def response_temperatures(self):
        """
        Returns a list of temperature grids (K) for each filter response.
        """
        return [r.temperature for r in self.responses]

    @property
    def response_values(self):
        """
        Returns a list of response values (DN cm^5 / pix / s) for each filter.
        """
        return [r.response for r in self.responses]

    @property
    def min_T(self):
        """
        Lower bound of log10 temperature grid.
        """
        return self._min_T

    @property
    def max_T(self):
        """
        Upper bound of log10 temperature grid.
        """
        return self._max_T

    @property
    def dT(self):
        """
        Bin width of log10 temperature grid.
        """
        return self._dT

    @property
    def min_error(self):
        """
        Minimum error applied to DN/s when intensity error is not provided.
        """
        return self._min_error

    @property
    def relative_error(self):
        """
        Relative error (%) used to scale intensity if error is not provided.
        """
        return self._relative_error

    @property
    def intensity_errors(self) -> u.Quantity:
        """
        Returns the intensity uncertainties, either user-provided or model-based.

        If not provided, errors are estimated using:
            max(relative_error * observed_intensity, min_error)

        For details, see:
        https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro

        Returns
        -------
        `~astropy.units.Quantity`
            Intensity errors in DN/s for each filter.
        """
        if self._intensity_errors is not None:
            return self._intensity_errors * (u.DN / u.s)

        warnings.warn(
            "No intensity_errors provided. Using default model: "
            f"max(relative_error * observed_intensity, min_error)\n"
            f"=> relative_error = {self.relative_error}, min_error = {self.min_error} DN/s\n"
            "See: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro",
            UserWarning,
        )

        estimated = np.maximum(
            self.relative_error * self._observed_intensities,
            self.min_error,
        )
        return estimated * (u.DN / u.s)

    @property
    def monte_carlo_runs(self) -> int:
        """
        Number of Monte Carlo runs to perform (0 = disabled).
        """
        return self._monte_carlo_runs

    @property
    def solv_factor(self):
        """
        Normalization factor used during DEM fitting to stabilize the solver.
        Default is 1e21.
        """
        return self._solv_factor

    @property
    def max_iterations(self):
        """
        Maximum number of iterations used in the least-squares DEM solver.
        """
        return self._max_iterations

    def create_logT_grid(self):
        """
        Build the DEM temperature grid *exactly* from min to max in steps of dT.
        """
        n_bins = int(round((self._max_T - self._min_T) / self._dT)) + 1
        self.logT = np.linspace(self._min_T, self._max_T, n_bins)
        self.T = (10**self.logT) * u.K
        self.dlogT = self._dT  # dimensionless - convenience-27
        self.dlnT = np.log(10.0) * self.dlogT  # needed for IDL-style integrals
        
    def _dem_per_log10T(self, dem_per_K):
        """Convert DEM per K → DEM per log10 T (cm^-5)."""

        return (np.log(10.0) * self.T) * dem_per_K


    def _interpolate_responses_to_grid(
        self,
    ):  # This mirrors what xrt_dem_iter_estim.pro does.
        """
        Interpolates each filter's temperature response onto the DEM temperature grid (self.logT).
        """
        self.interpolated_responses = []

        for i, (T_orig, R_orig, fname) in enumerate(
            zip(self.response_temperatures, self.response_values, self.filter_names)
        ):
            logT_orig = np.log10(T_orig.to_value(u.K))
            response_vals = R_orig.to_value(u.DN * u.cm**5 / (u.pix * u.s))

            interp_func = interp1d(
                logT_orig,
                response_vals,
                kind="linear",
                bounds_error=False,
                fill_value=0.0,
            )  # kind = 'cubic' )  kind="linear",

            R_interp = interp_func(self.logT)
            self.interpolated_responses.append(R_interp)

    @property
    def response_matrix(self):
        """
        Returns the response matrix after interpolation.

        Shape: (n_filters, n_temperatures)
        """
        if not hasattr(self, "_response_matrix"):
            raise AttributeError(
                "Response matrix has not been built yet. Call _build_response_matrix()."
            )
        return self._response_matrix

    def _build_response_matrix(self):
        """
        Builds the response matrix from interpolated responses.

        Sets:
        -------
        self.response_matrix : ndarray
            2D array of shape (n_filters, n_temperatures)
            Stack your self.interpolated_responses into a 2D NumPy array

        Personal notes: The response matrix is a 2D array that relates temperature to observed intensity
        For numerical DEM:
            -You approximate the integral as a matrix multiplication
            -Each filter contributes one equation (row)
            -Each temperature bin contributes one unknown (column)
            - Intergal DEM(T) * R(T) dT = sum[DEM_i * R_i * dT]
        """
        if not hasattr(self, "interpolated_responses"):
            raise RuntimeError(
                "Call _interpolate_responses_to_grid() before building the response matrix."
            )

        # self._response_matrix = np.vstack(self.interpolated_responses) # matrix
        self._response_matrix = np.vstack(self.interpolated_responses).astype(
            float
        )  # matrix

        print(
            f"Built response matrix: shape = {self._response_matrix.shape} (filters * logT bins)"
        )

    def _estimate_initial_dem(
        self, smooth=False, logscale=False, plot=True
    ):  # mirrors xrt_dem_iter_estim.pro
        """
        Estimates an initial DEM by inverting I ≈ R @ DEM.
        An initial estimate of how much plasma is emitting at each temperature

        Sets
        -----
        self.initial_dem : np.ndarray
            First-guess DEM estimate across logT grid (length = n_temperatures)

        Parameters
        ----------
        smooth : bool
            If True, applies mild Gaussian smoothing (future option).
        logscale : bool
            If True, computes log10(DEM) and exponentiates to suppress spikes.
        plot : bool
            If True, shows a diagnostic plot of the initial DEM.
        """
        print(self.response_temperatures, self.response_values, self.filter_names)

        # Define inputs - xrt_dem_iter_estim.pro
        if not hasattr(self, "response_matrix"):
            raise RuntimeError("Run _build_response_matrix() before estimating DEM.")

        I_obs = np.asarray(self._observed_intensities, dtype=float)
        R = self.response_matrix
        n_filters, n_temps = R.shape

        print(
            f"Estimating DEM from {n_filters} intensities across {n_temps} temperature bins..."
        )

        with np.errstate(divide="ignore", invalid="ignore"):
            # First estimate DEM per-logT (cm^-5)
            dem_logT = np.zeros(n_temps) #27
            for i in range(n_filters):
                row = R[i, :]
                ratio = np.where(row > 1e-30, I_obs[i] / row, 0.0)  # cm^-5
                dem_logT += ratio
                
            dem_logT /= n_filters
            
        # DO NOT divide by self._dT here.
        # Optional smoothing in per-logT space:
        if smooth:
            from scipy.ndimage import gaussian_filter1d
            dem_logT = gaussian_filter1d(dem_logT, sigma=1.0)

        # Store canonical DEM (per-logT) for fitting/forward model
        self.initial_dem_logT = dem_logT * (u.cm**-5) ##################################
        #self.initial_dem_logT = (np.log(10.0) * self.T) * self.dem_initial * u.cm**-5  # (N,), per-log10T #28

        # For plotting PER-K if you want that axis:
        dem_perK = (dem_logT / (np.log(10.0) * self.T.to_value(u.K)))  # cm^-5 K^-1
        self.initial_dem = dem_perK * (u.cm**-5 / u.K)
        
        #self.dem_initial = dem_per_K_on_grid  # (N,), per-K (cm^-5 K^-1), np.ndarray or Quantity    
        
    
            #estimates = np.zeros(n_temps)
            # for i in range(n_filters):
            #     row = R[i, :]
            #     ratio = np.where(row > 1e-30, I_obs[i] / row, 0.0)
            #     estimates += ratio

            # estimates /= n_filters
            # estimates /= self._dT  # Convert to per-logT-bin definition
            # assert estimates.shape[0] == len(self.logT)

            # if logscale:
            #     # Suppress large dynamic range and spikes
            #     estimates = 10 ** np.log10(estimates + 1e-30)

            # if smooth:
            #     from scipy.ndimage import gaussian_filter1d

            #     estimates = gaussian_filter1d(estimates, sigma=1.0)

        # estimates =dem_logT
        # # Apply units
        # self.initial_dem = estimates * (u.cm**-5 / u.K)
        # print("  Max:", np.max(estimates))
        # print("  Min:", np.min(estimates))
        # print("  dT (logT bin size):", self._dT)

        print("  Max DEM_per_logT:", np.max(dem_logT))
        print("  Min DEM_per_logT:", np.min(dem_logT))
        print("  dlog10T:", self.dlogT)

        print("Initial DEM estimate complete")
        print(f"Peak DEM_per_logT: {self.initial_dem_logT.max():.2e}")
        print(f" Mean DEM_per_logT: {self.initial_dem_logT.mean():.2e}")

        # Diagnostics

        print(f"I_obs: {I_obs}")  # Observed intensities
        print(f"R (response matrix): {R.shape}")
        print(f"Sum of response rows: {[np.sum(R[i]) for i in range(R.shape[0])]}")
        print(f"dT: {self._dT}")
        print("[DEBUG] DEM before dT division:")

        print("[DEBUG] Response row sums:")
        for i, row in enumerate(R):
            print(
                f"  {self.filter_names[i]}: sum={np.sum(row):.2e}, max={np.max(row):.2e}"
            )

        print(f"[DEBUG] dT: {self._dT:.3f}")

        # Plotting
        if plot:
            import matplotlib.pyplot as plt

            plt.figure(figsize=(8, 4))
            ylabel = "DEM [cm⁻⁵ K⁻¹]"

            # Custom label with filters and date
            filters_str = ", ".join(self.observed_channel)
            label_str = f"Initial DEM\n{filters_str}\n"  # {self.date_obs}"

            if logscale:
                plt.plot(self.logT, np.log10(self.initial_dem.value), drawstyle="steps-mid",
                        label=label_str, color="purple")
                ylabel = r"log$_{10}$ DEM [cm$^{-5}$ K$^{-1}$]"
            else:
                plt.plot(self.logT, self.initial_dem.value, drawstyle="steps-mid",
                        label=label_str, color="purple")
                plt.yscale("log")

            # if logscale:
            #     log_dem_vals = np.log10(self.initial_dem.value + 1e-30)
            #     plt.plot(
            #         self.logT,
            #         log_dem_vals,
            #         drawstyle="steps-mid",
            #         label=label_str,
            #         color="purple",
            #     )
            #     ylabel = "log₁₀ DEM [cm⁻⁵ K⁻¹]"
            # else:
            #     plt.plot(
            #         self.logT,
            #         self.initial_dem.value,
            #         drawstyle="steps-mid",
            #         label=label_str,
            #         color="purple",
            #     )
            #     plt.yscale("log")

            plt.xlabel("log₁₀ T [K]")
            plt.ylabel(ylabel)
            plt.title("Initial DEM Estimate")
            plt.grid(True)
            plt.legend(loc="upper right", fontsize=8)
            plt.tight_layout()
            plt.show()

        print("Initial DEM estimate complete")

    # STEP 1 - Each temperature bin gets its own parameter, initialized with your initial DEM estimate
    def _build_lmfit_parameters(self):
        """
        Initializes lmfit Parameters from the initial DEM guess.

        Sets:
        -------
        self.lmfit_params : lmfit.Parameters
            Each temperature bin gets a parameter (free by default).
        """
        # if not hasattr(self, "initial_dem"):
        #     raise RuntimeError(
        #         "Call _estimate_initial_dem() before building parameters."
        #     )

        # params = Parameters()

        # for i, val in enumerate(self.initial_dem):
        #     # You could add bounds here if needed (e.g., min=0)
        #     params.add(f"dem_{i}", value=val, min=0)
        
        #27
        # for i, val in enumerate(self.initial_dem):
        #     # Convert to float if it's a Quantity
        #     if hasattr(val, "unit"):
        #         val = val.to_value()  # default: returns value in current unit
        #     params.add(f"dem_{i}", value=val, min=0)

        # self.lmfit_params = params
        # print(f"Built {len(params)} lmfit parameters for DEM fit")



        if not hasattr(self, "initial_dem_logT"):
            raise RuntimeError("Call _estimate_initial_dem() before building parameters.")
        params = Parameters()
        for i, val in enumerate(self.initial_dem_logT.to_value(u.cm**-5)):
            params.add(f"dem_{i}", value=float(val), min=0.0)
        self.lmfit_params = params
        print(f"Built {len(params)} lmfit parameters for DEM fit")


    # STEP 2: Build the residual function
    # This function computes how far off your DEM model’s predicted intensities are from your observed ones, normalized by the uncertainty.
    def _residuals(self, params):
        """
        Computes the residuals between modeled and observed intensities.

        Parameters
        ----------
        params : lmfit.Parameters
            DEM values at each temperature bin.

        Returns
        -------
        np.ndarray
            Residuals = (I_model - I_obs) / sigma
        """
        # # 1. Get DEM vector from lmfit Parameters
        # dem_vector = np.array([params[f"dem_{i}"].value for i in range(len(self.logT))])

        # DEM per-logT (cm^-5)
        dem_logT = np.array([params[f"dem_{i}"].value for i in range(len(self.logT))])


        # 2. Compute modeled intensities: I_model = R · DEM
        #I_model = self.response_matrix @ dem_vector
        I_model = self.response_matrix @ (dem_logT * self.dlogT) #27

        # 3. Determine observational errors (user-provided or fallback)
        if self._intensity_errors is not None:
            errors = np.array(self._intensity_errors)
        else:
            errors = np.maximum(
                self.min_error, self.relative_error * self._observed_intensities
            )

        # 4. Return normalized residuals
        residuals = (I_model - self._observed_intensities) / errors
        print(
            "[•] Residuals stats → mean: {:.2e}, std: {:.2e}".format(
                np.mean(residuals), np.std(residuals)
            )
        )
        return residuals

    def fit_dem(self):
        """
        Runs the DEM fitting using lmfit's least-squares minimization.

        Sets:
        -------
        self.fitted_dem : np.ndarray
            Best-fit DEM solution (length = n_temps)
        self.result : lmfit.MinimizerResult
            Full fit result object from lmfit
        """
        # if not hasattr(self, "lmfit_params"):
        #     self._build_lmfit_parameters()

        # if not hasattr(self, "lmfit_params"):
        #     raise RuntimeError("Call _build_lmfit_parameters() before fitting.")

        #27
        if not hasattr(self, "initial_dem_logT"):
            raise RuntimeError("Call _estimate_initial_dem() first.")
        params = Parameters()
        for i, val in enumerate(self.initial_dem_logT.to_value(u.cm**-5)):
            params.add(f"dem_{i}", value=float(val), min=0.0)
        self.lmfit_params = params

        print("Starting DEM optimization..")
        result = minimize(
            self._residuals,
            self.lmfit_params,
            method="least_squares",
            max_nfev=self.max_iterations,
        )

        self.result = result
        
        dem_best_logT = np.array([self.result.params[f"dem_{i}"].value for i in range(len(self.logT))])

        self.fitted_dem_logT = dem_best_logT * (u.cm**-5) #28
        self.fitted_dem = (dem_best_logT / (np.log(10.0) * self.T.to_value(u.K))) * (u.cm**-5 / u.K)
        
        #self.dem_fit = dem_per_K_on_grid  # (N,), per-K (cm^-5 K^-1)
        #self.fitted_dem_logT = (np.log(10.0) * self.T) * self.dem_fit * u.cm**-5  # (N,), per-log10T

        if not result.success:
            print("[⚠️] DEM fit did not fully converge:")
            print("  →", result.message)

        # self.fitted_dem = np.array([result.params[f"dem_{i}"].value for i in range(len(self.logT))])

        print(
            f"[✓] DEM fit complete — reduced chi-squared: {result.chisqr / len(self._observed_intensities):.2f}"
        )

        print("[✓] DEM fit complete")
        print(
            f"  → Reduced chi-squared: {result.chisqr / len(self._observed_intensities):.2f}"
        )
        print(f"  → Total iterations: {result.nfev}")

        return result

    def print_residual_diagnostics(self, params):

        dem_logT = np.array([params[f"dem_{i}"].value for i in range(len(self.logT))])
        I_model = self.response_matrix @ (dem_logT * self.dlogT)
        residuals = (I_model - self._observed_intensities) / self._intensity_errors
        
        print("Observed Intensities:", self._observed_intensities)
        print("Modeled Intensities:", I_model)
        print("Errors:", self._intensity_errors)
        print("Residuals:", residuals)
        print(
            f"[•] Residuals stats → mean: {residuals.mean():.2e}, std: {residuals.std():.2e}"
        )


    def plot_dem_fit(
        self,
        logscale: bool = True,
        scale: str = "per_K",        # "per_K" or "per_log10T"
        show_initial: bool = True,
        ax=None,
        title: str | None = None,
    ):
        """
        Plot the fitted DEM (and optional initial DEM) using a consistent scale.

        logscale=True  -> semilogy of linear values
        logscale=False -> linear plot of log10(values)

        scale="per_K"       -> DEM [cm^-5 K^-1]
        scale="per_log10T"  -> phi = DEM * T * ln(10) [cm^-5]
        """
        import numpy as np
        import matplotlib.pyplot as plt

        if ax is None:
            fig, ax = plt.subplots()

        # Grid
        logT = self.logT
        T = 10.0 ** logT

        # --- find fitted DEM (prefer canonical names; fall back to legacy) ---
        linear_candidates = ["dem_fit", "dem", "fitted_dem", "dem_solved", "dem_solution"]
        log_candidates    = ["logdem_fit", "logdem", "fitted_logdem", "logdem_solved"]

        dem_fit_perK = None
        for name in linear_candidates:
            if hasattr(self, name) and getattr(self, name) is not None:
                dem_fit_perK = np.asarray(getattr(self, name), dtype=float)
                break
        if dem_fit_perK is None:
            for name in log_candidates:
                if hasattr(self, name) and getattr(self, name) is not None:
                    dem_fit_perK = 10.0 ** np.asarray(getattr(self, name), dtype=float)
                    break
        if dem_fit_perK is None:
            raise RuntimeError(
                "No fitted DEM found. Expected one of "
                f"{linear_candidates + log_candidates} to exist on the object."
            )

        # --- initial DEM (optional) ---
        dem_init_perK = None
        if show_initial:
            init_linear_candidates = ["initial_dem", "dem_initial"]
            init_log_candidates    = ["initial_logdem", "logdem_initial"]
            for name in init_linear_candidates:
                if hasattr(self, name) and getattr(self, name) is not None:
                    dem_init_perK = np.asarray(getattr(self, name), dtype=float)
                    break
            if dem_init_perK is None:
                for name in init_log_candidates:
                    if hasattr(self, name) and getattr(self, name) is not None:
                        dem_init_perK = 10.0 ** np.asarray(getattr(self, name), dtype=float)
                        break

        # --- choose scientific scale ---
        if scale == "per_K":
            y_fit_linear  = np.clip(dem_fit_perK, np.finfo(float).tiny, None)
            y_init_linear = None if dem_init_perK is None else np.clip(dem_init_perK, np.finfo(float).tiny, None)
            y_label_lin   = r"DEM per K  [cm$^{-5}$ K$^{-1}$]"
            y_label_log10 = r"$\log_{10}$ DEM per K  [cm$^{-5}$ K$^{-1}$]"
        elif scale == "per_log10T":
            y_fit_linear  = np.clip(dem_fit_perK * T * np.log(10.0), np.finfo(float).tiny, None)
            y_init_linear = None if dem_init_perK is None else np.clip(dem_init_perK * T * np.log(10.0), np.finfo(float).tiny, None)
            y_label_lin   = r"DEM per $\log_{10}T$  [cm$^{-5}$]"
            y_label_log10 = r"$\log_{10}$ DEM per $\log_{10}T$  [cm$^{-5}$]"
        else:
            raise ValueError("scale must be 'per_K' or 'per_log10T'")

        # --- plot without double-logging ---
        if logscale:
            ax.semilogy(logT, y_fit_linear, label="Fitted DEM")
            if y_init_linear is not None:
                ax.semilogy(logT, y_init_linear, "--", alpha=0.7, label="Initial DEM")
            ax.set_ylabel(y_label_lin)
        else:
            ax.plot(logT, np.log10(y_fit_linear), label="Fitted DEM")
            if y_init_linear is not None:
                ax.plot(logT, np.log10(y_init_linear), "--", alpha=0.7, label="Initial DEM")
            ax.set_ylabel(y_label_log10)

        ax.set_xlabel(r"$\log_{10} T$ [K]")
        ax.set_title(title or ("Initial vs Fitted DEM" if show_initial else "Fitted DEM"))
        ax.legend()
        ax.grid(True, alpha=0.3)


    # def plot_dem_fit(self, logscale=True):
    #     import matplotlib.pyplot as plt
    #     # if not hasattr(self, "initial_dem_logT"):
    #     #     raise RuntimeError("Initial DEM not computed. Run _estimate_initial_dem() first.")
    #     # if not hasattr(self, "result"):
    #     #     raise RuntimeError("DEM fit result not available. Run fit_dem() first.")


    #     # --- grid sanity ---
    #     if not hasattr(self, "logT"):
    #         raise RuntimeError("Temperature grid missing. Run create_logT_grid() first.")
    #     # ensure we have linear T for conversion if needed
    #     if not hasattr(self, "T"):
    #         self.T = np.power(10.0, self.logT)
        
    #     ln10T = np.log(10.0) * self.T
        
    #     # Choose per-logT plotting (recommended for XRT)
    #     initial_vals = self.initial_dem_logT.to_value(u.cm**-5)
    #     best_vals    = self.fitted_dem_logT.to_value(u.cm**-5)
    #     ylabel = r"DEM per $\log_{10}T$ [cm$^{-5}$]"

    #     plt.figure(figsize=(10, 5))
    #     plt.plot(self.logT, initial_vals, drawstyle="steps-mid", label="Initial DEM", linestyle="--", color="gray")
    #     plt.plot(self.logT, best_vals,    drawstyle="steps-mid", label="Fitted DEM",  color="blue")

    #     plt.xlabel("log₁₀ T [K]")
    #     plt.title("Initial vs Fitted DEM")
    #     plt.legend()
    #     plt.grid(True)
    #     if logscale:
    #         plt.yscale("log")
    #     plt.ylabel(ylabel)
    #     plt.tight_layout()
    #     plt.show()

    #     print(f"[Plot] Peak Initial DEM (per-logT): {np.max(initial_vals):.2e}")
    #     print(f"[Plot] Peak Fitted  DEM (per-logT): {np.max(best_vals):.2e}")


    # def plot_dem_fit(self, logscale=True):
    #     """
    #     Plots the initial and fitted DEM on the same logT grid.

    #     Parameters
    #     ----------
    #     logscale : bool
    #         If True, uses a logarithmic y-axis.
    #     """
    #     import matplotlib.pyplot as plt

    #     if not hasattr(self, "initial_dem"):
    #         raise RuntimeError(
    #             "Initial DEM not computed. Run _estimate_initial_dem() first."
    #         )
    #     if not hasattr(self, "result"):
    #         raise RuntimeError("DEM fit result not available. Run fit_dem() first.")

    #     # Extract best-fit DEM from lmfit result
    #     # best_fit_vals = np.array(
    #     #     [self.result.params[f"dem_{i}"].value for i in range(len(self.logT))]
    #     # )
    #     # initial_dem_vals = (
    #     #     self.initial_dem.value
    #     #     if hasattr(self.initial_dem, "value")
    #     #     else self.initial_dem
    #     # )
    #     # log_initial_dem_vals = np.log10(self.initial_dem.value) if hasattr( np.log10(self.initial_dem), "value") else np.log10(self.initial_dem)

    #     initial_vals = self.initial_dem_logT.to_value(u.cm**-5)
    #     best_vals    = self.fitted_dem_logT.to_value(u.cm**-5)

    #     plt.figure(figsize=(10, 5))
    #     plt.plot(
    #         self.logT,
    #         initial_dem_vals,
    #         drawstyle="steps-mid",
    #         label="Initial DEM",
    #         linestyle="--",
    #         color="gray",
    #     )
    #     plt.plot(
    #         self.logT,
    #         best_fit_vals,
    #         drawstyle="steps-mid",
    #         label="Fitted DEM",
    #         color="blue",
    #     )

    #     plt.xlabel("log₁₀ T [K]")
    #     # plt.ylabel("DEM [cm⁻⁵ K⁻¹]")
    #     ylabel = r"DEM per $\log_{10}T$ [cm$^{-5}$]"
    #     plt.title("Initial vs Fitted DEM")
    #     plt.legend()
    #     plt.grid(True)
    #     if logscale:
    #         plt.yscale("log")
    #         plt.ylabel(r"DEM [cm$^{-5}$ K$^{-1}$] (log-scaled)")
    #     else:
    #         plt.ylabel(r"DEM [cm$^{-5}$ K$^{-1}$]")
    #     plt.tight_layout()
    #     plt.show()

    #     print(f"[Plot] Peak Initial DEM: {np.max(initial_dem_vals):.2e}")
    #     print(f"[Plot] Peak Fitted DEM: {np.max(best_fit_vals):.2e}")

    def summary(self):
        print("XRTpy DEM Iterative Setup Summary")
        print("-" * 40)
        print(f" Filters:           {self.filter_names}")
        print(f" Obs Intensities:   {self.observed_intensities}")
        print(f" Number of observations (Nobs): {len(self._observed_intensities)}")
        print(f" Solver Normalization Factor: {self.solv_factor:.1e}")
        print(
            f" Monte Carlo runs:  {self.monte_carlo_runs if self.monte_carlo_runs > 0 else 'None'}"
        )
        print(f" Max Iterations:    {self.max_iterations}")
        print(f" Intensity Errors:  {self.intensity_errors}")
        print(f" Temp Grid:         logT {self.min_T} to {self.max_T} (step {self.dT})")
        print(f" Temp bins:         {len(self.logT)}")
        print(
            f" Error model used:  {'User-provided' if self._intensity_errors is not None else 'Auto (obs * 0.03, min=2 DN/s)'}"
        )
        if self._intensity_errors is None:
            print(
                "For more info: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro"
            )
        print("-" * 40)

__all__ = [
    "XRTDEMIterative",
]

import warnings

import astropy.units as u
import numpy as np

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
        solv_factor=1e21
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
        self.observed_channel = validate_and_format_filters(observed_channel)

        if observed_channel is None or len(observed_channel) == 0:
            raise ValueError("`observed_channel` is required and cannot be empty.")

        # Store intensity and error arrays
        self._observed_intensities = np.asarray(observed_intensities, dtype=float)

        if observed_intensities is None or len(observed_intensities) == 0:
            raise ValueError("`observed_intensities` is required and cannot be empty.")

        # Errors
        if intensity_errors is not None:
            self._intensity_errors = np.asarray(intensity_errors, dtype=float)
            if self._intensity_errors.shape != self.observed_intensities.shape:
                raise ValueError(
                    "Length of intensity_errors must match observed_intensities."
                )
        else:
            self._intensity_errors = None  # Will be computed later

        # Store temperature grid parameters
        self._min_T = float(min_T)
        self._max_T = float(max_T)
        self._dT = float(dT)

        # Validate Monte Carlo setting
        if isinstance(monte_carlo_runs, bool):
            raise ValueError(
                "monte_carlo_runs must be a non-negative whole number, not a boolean."
            )
        elif isinstance(monte_carlo_runs, (int, np.integer)) or isinstance(monte_carlo_runs, float) and monte_carlo_runs.is_integer():
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


    def _interpolate_responses_to_grid(self): #This mirrors what xrt_dem_iter_estim.pro does.
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
            )

            R_interp = interp_func(self.logT)
            self.interpolated_responses.append(R_interp)




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
            raise RuntimeError("Call _interpolate_responses_to_grid() before building the response matrix.")
        
        #self._response_matrix = np.vstack(self.interpolated_responses) # matrix
        self._response_matrix = np.vstack(self.interpolated_responses).astype(float) # matrix
        
        print(f"Built response matrix: shape = {self._response_matrix.shape} (filters * logT bins)")



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

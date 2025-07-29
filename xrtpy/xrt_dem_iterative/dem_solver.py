__all__ = [
    "XRTDEMIterative",
]

import astropy.time
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters , minimize
from scipy.interpolate import interp1d, CubicSpline
import warnings
from numpy import trapezoid #np.trapz Deprecation

from xrtpy.util.filters import solve_filter_name, validate_and_format_filters
from xrtpy.util.time import epoch



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


        #Errors
        if intensity_errors is not None:
            self._intensity_errors = np.asarray(intensity_errors, dtype=float)
            if self._intensity_errors.shape != self.observed_intensities.shape:
                raise ValueError("Length of intensity_errors must match observed_intensities.")
        else:
            self._intensity_errors = None  # Will be computed later
            
        
        # Store temperature grid parameters
        self._min_T = float(min_T)
        self._max_T = float(max_T)
        self._dT = float(dT)
        
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
                raise ValueError("Length of intensity_errors must match observed_intensities.")
        else:
            self._intensity_errors = None
    

    def __repr__(self):
        return f"<XRTDEMIterative(filters={self.filter_names}, logT={self.min_T}–{self.max_T}, dT={self.dT})>"
        
    # @property  #Removed if not used
    # def name(self) -> str:
    #     """
    #     The XRT filter channel name, standardized (e.g. "Al-mesh").
    #     """
    #     return self._name

    @property
    def observed_intensities(self) -> u.Quantity: #Add method to account for known values not worth observed_intensities
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
        """Lower bound of log10 temperature grid."""
        return self._min_T

    @property
    def max_T(self):
        """Upper bound of log10 temperature grid."""
        return self._max_T

    @property
    def dT(self):
        """Bin width of log10 temperature grid."""
        return self._dT
    
    @property
    def min_error(self):
        """Minimum error applied to DN/s when intensity error is not provided."""
        return self._min_error

    @property
    def relative_error(self):
        """Relative error (%) used to scale intensity if error is not provided."""
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
            UserWarning)

        estimated = np.maximum(
            self.relative_error * self._observed_intensities,
            self.min_error,
        )
        return estimated * (u.DN / u.s)


    def summary(self):
        print("XRTpy DEM Iterative Setup Summary")
        print("-" * 40)
        print(f" Filters:           {self.filter_names}")
        print(f" Obs Intensities:   {self.observed_intensities}")
        print(f" Number of observations (Nobs): {len(self._observed_intensities)}")
        print(f" Intensity Errors:  {self.intensity_errors}")
        print(f" Temp Grid:         logT {self.min_T} to {self.max_T} (step {self.dT})")
        print(f" Temp bins:         {len(self.logT)}")
        print(f" Error model used:  {'User-provided' if self._intensity_errors is not None else 'Auto (obs * 0.03, min=2 DN/s)'}")
        if self._intensity_errors is None:
            print("For more info: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro")
        print("-" * 40)

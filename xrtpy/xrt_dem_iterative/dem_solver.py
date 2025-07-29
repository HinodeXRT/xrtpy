__all__ = [
    "XRTDEMIterative",
]

import astropy.time
import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np
from lmfit import Parameters , minimize
from scipy.interpolate import interp1d, CubicSpline

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
        # Validate and store filter names
        self.observed_channel = validate_and_format_filters(observed_channel)

        # Store intensity and error arrays
        self._observed_intensities = np.asarray(observed_intensities, dtype=float)

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

        print(
            "\n[INFO] No intensity_errors provided. "
            "Using default model: max(relative_error * observed_intensity, min_error)\n"
            "       => relative_error = {:.2f}, min_error = {:.1f} DN/s\n"
            "       => For details: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro\n".format(
                self.relative_error, self.min_error
            )
        )

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
        print(f" Intensity Errors:  {self.intensity_errors}")
        print(f" Temp Grid:         logT {self.min_T} to {self.max_T} (step {self.dT})")
        print(f" Temp bins:         {len(self.logT)}")
        print(f" Error model used:  {'User-provided' if self._intensity_errors is not None else 'Auto (obs * 0.03, min=2 DN/s)'}")
        if self._intensity_errors is None:
            print("For more info: https://hesperia.gsfc.nasa.gov/ssw/hinode/xrt/idl/util/xrt_dem_iterative2.pro")
        print("-" * 40)

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
        self.observed_intensities = np.asarray(observed_intensities, dtype=float)
        if intensity_errors is not None:
            self.intensity_errors = np.asarray(intensity_errors, dtype=float)
            if self.intensity_errors.shape != self.observed_intensities.shape:
                raise ValueError("Length of intensity_errors must match observed_intensities.")
        else:
            self.intensity_errors = None  # Will be computed later

        # Store temperature response objects
        self.responses = temperature_responses

        # Store temperature grid parameters
        self.min_T = min_T
        self.max_T = max_T
        self.dT = dT
        self.logT = np.arange(min_T, max_T + dT / 2, dT)

        # Store error model parameters
        self.min_error = min_error
        self.relative_error = relative_error
        
    
    @property  #Removed if not used
    def name(self) -> str:
        """
        The XRT filter channel name, standardized (e.g. "Al-mesh").
        """
        return self._name

    @property
    def observed_intensities(self) -> u.Quantity:
        """
        Observed intensities with physical units.
        Returns
        -------
        `~astropy.units.Quantity`
            Intensities in DN/s for each filter channel.
        """
        return self._observed_intensities * (u.DN / u.s)
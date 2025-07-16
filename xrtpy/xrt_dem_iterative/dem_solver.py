__all__ = [
#    "",
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



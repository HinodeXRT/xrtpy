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
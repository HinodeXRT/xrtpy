__all__ = [
    "ComputeDEMStatistics",
]

import numpy as np


class ComputeDEMStatistics:
    """
    Diagnostic class for computing residual statistics from a fitted DEM solution.

    This class provides utilities to:
        Compute chi-squared and reduced chi-squared between observed and modeled intensities.
        Print residual diagnostics per filter (similar to IDL's xrt_iter_demstat.pro).

    Methods
    -------
    compute_chi_squared()
        Compute total and reduced chi-squared for DEM fit.

    print_residuals()
        Print modeled vs. observed intensities and residuals (normalized by error).
    """
    
    def __init__(self, dem_solver):
        self.dem_solver = dem_solver


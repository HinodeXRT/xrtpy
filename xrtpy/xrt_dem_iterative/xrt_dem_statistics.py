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

    def compute_chi_squared(self):
        """
        Compute chi-squared and reduced chi-squared between observed and modeled intensities.

        Returns
        -------
        chi2 : float
            Total chi-squared value.
        chi2_red : float or None
            Reduced chi-squared (chi2 / dof), or None if dof <= 0.
        """
        if not hasattr(self.dem_solver, "fitted_dem"):
            raise RuntimeError("Must run fit_dem() before computing chi-squared.")

        I_model = self.dem_solver.response_matrix @ self.dem_solver.fitted_dem
        I_obs = self.dem_solver._observed_intensities

        abs_error = np.maximum(
            self.dem_solver.min_error,
            self.dem_solver.relative_error * I_obs
        )

        chi2 = np.sum(((I_model - I_obs) / abs_error) ** 2)
        dof = len(I_obs) - len(self.dem_solver.fitted_dem)
        chi2_red = chi2 / dof if dof > 0 else None

        return chi2, chi2_red

    def print_residuals(self):
        """
        Print residuals and modeled vs. observed intensities (IDL-style diagnostics).
        """
        I_model = self.dem_solver.response_matrix @ self.dem_solver.fitted_dem
        I_obs = self.dem_solver._observed_intensities
        abs_error = np.maximum(
            self.dem_solver.min_error,
            self.dem_solver.relative_error * I_obs
        )

        residuals = (I_model - I_obs) / abs_error

        print("\n[DEM RESIDUALS PER FILTER]")
        print("---------------------------")
        for i, name in enumerate(self.dem_solver.filter_names):
            print(
                f"{name:<20}  Obs: {I_obs[i]:.2f}  Model: {I_model[i]:.2f}  "
                f"Error: {abs_error[i]:.2f}  Residual: {residuals[i]:+.2f}"
            )
        print(f"\nMean residual: {np.mean(residuals):+.2f}")
        print(f"Std  residual: {np.std(residuals):.2f}")

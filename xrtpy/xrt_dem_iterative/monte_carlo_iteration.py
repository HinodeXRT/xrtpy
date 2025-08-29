__all__ = [
    "MonteCarloIteration",
]

import numpy as np


class MonteCarloIteration:
    
    def __init__(self, dem_solver):
        """
        Parameters
        ----------
        dem_solver : XRTDEMIterative
            A fitted DEM object with observed intensities, errors, and temperature grid.
        """
        self.dem_solver = dem_solver

        if not hasattr(dem_solver, "logT"):
            raise RuntimeError("DEM solver must have a defined temperature grid.")
        if not hasattr(dem_solver, "intensity_errors"):
            raise RuntimeError("DEM solver must define intensity errors.")

        self.n_bins = len(dem_solver.logT)
        self.n_filters = len(dem_solver.observed_intensities)

    def generate_mc_realizations(self, n_realizations=100, seed=None,reject_negative=True):
        """
        Generate randomized intensity realizations for Monte Carlo uncertainty estimation.

        Parameters
        ----------
        n_realizations : int
            Number of Monte Carlo runs to generate.
        seed : int or None
            Random seed for reproducibility.

        Sets
        ----
        self.mc_intensity_sets : np.ndarray
            Shape (n_realizations, n_filters), randomized intensities.
        """
        if seed is not None:
            np.random.seed(seed)

        # Compute error bars
        abs_error = np.maximum(
            self.min_error, self.relative_error * self._observed_intensities
        )

        # Draw random perturbations for each intensity
        self.mc_intensity_sets = np.random.normal(
            loc=self._observed_intensities,
            scale=abs_error,
            size=(n_realizations, len(self._observed_intensities)),
        )

    def run_mc_simulation(self, n_realizations=100, seed=None):
        """
        Run Monte Carlo simulations to estimate DEM uncertainties.

        Parameters
        ----------
        n_realizations : int
            Number of Monte Carlo realizations to run.
        seed : int or None
            Optional seed for reproducibility.

        Sets
        ----
        self.mc_dems : np.ndarray
            Shape (n_temps, n_realizations). Each column is a DEM realization.
        """
        if seed is not None:
            np.random.seed(seed)

        # Use user-provided or fallback error model
        if self._intensity_errors is not None:
            errors = np.array(
                self._intensity_errors, dtype=float
            )  # Covering given user error in pyton array
            # errors = self._intensity_errors
        else:
            errors = np.maximum(
                self.min_error, self.relative_error * self._observed_intensities
            )

        self.mc_intensity_sets = np.random.normal(
            loc=self._observed_intensities[:, None],  # shape (5, 1)
            scale=errors[:, None],  # shape (5, 1)
            size=(len(self._observed_intensities), n_realizations),  # shape (5, 20)
        )
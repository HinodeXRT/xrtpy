__all__ = [
    "Monte_Carlo_Iteration",
]

import numpy as np


class Monte_Carlo_Iteration:

    # def __init__( ):

    def generate_mc_realizations(self, n_realizations=100, seed=None):
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

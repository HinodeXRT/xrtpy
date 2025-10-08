__all__ = [
    "plot_dem_results",
    "plot_dem_uncertainty",
    "plot_idl_style",
    "plot_fit_residuals",
    "plot_dem_with_median_bins",
    "plot_iteration_stats",
]

import matplotlib.pyplot as plt
import numpy as np
import astropy.units as u


###############################################################################
# Plotting function 1
###############################################################################
def plot_dem_results(dem):
    """
    Plot the fitted DEM solution (with optional Monte Carlo uncertainty).

    Parameters
    ----------
    dem : XRTDEMIterative
        Solver object. If not yet solved, .solve() will be called.
    """
    if not hasattr(dem, "dem"):
        dem.solve()

    logT = dem.logT
    best_fit = dem.dem
    dem_err = getattr(dem, "dem_uncertainty", None)

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.step(
        logT,
        np.log10(best_fit + 1e-40),
        where="mid",
        color="blue",
        linewidth=2,
        label="Best-fit DEM",
    )

    if dem_err is not None:
        upper = np.log10(best_fit + dem_err + 1e-40)
        lower = np.log10(np.clip(best_fit - dem_err, 1e-40, None))
        ax.fill_between(
            logT, lower, upper, step="mid", color="blue", alpha=0.2, label="+/-1σ"
        )

    ax.set_xlabel("log10 T [K]")
    ax.set_ylabel("log10 DEM [cm$^{-5}$ K$^{-1}$]")
    ax.set_xlim(logT.min(), logT.max())
    ax.set_ylim(
        np.floor(np.log10(best_fit.min() + 1e-40)),
        np.ceil(np.log10(best_fit.max() + 1e-40)),
    )
    ax.set_title("DEM Solution")
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


###############################################################################
# Plotting function 2
###############################################################################
def plot_dem_uncertainty(dem):
    """
    Plot DEM with Monte Carlo uncertainty band.
    """
    if not hasattr(dem, "dem"):
        dem.solve()
    if not hasattr(dem, "dem_uncertainty"):
        raise AttributeError("No DEM uncertainty found. Run with monte_carlo_runs > 0.")

    logT = dem.logT
    best_fit = dem.dem
    dem_err = dem.dem_uncertainty

    fig, ax = plt.subplots(figsize=(8, 6))

    ax.step(
        logT,
        np.log10(best_fit + 1e-40),
        where="mid",
        color="blue",
        linewidth=2,
        label="Best-fit DEM",
    )

    upper = np.log10(best_fit + dem_err + 1e-40)
    lower = np.log10(np.clip(best_fit - dem_err, 1e-40, None))
    ax.fill_between(
        logT, lower, upper, step="mid", color="blue", alpha=0.2, label="+/-1σ"
    )

    ax.set_xlabel("log10 T [K]")
    ax.set_ylabel("log10 DEM [cm$^{-5}$ K$^{-1}$]")
    ax.set_xlim(logT.min(), logT.max())
    ax.set_ylim(
        np.floor(np.log10(best_fit.min() + 1e-40)),
        np.ceil(np.log10((best_fit + dem_err).max() + 1e-40)),
    )
    ax.set_title("DEM with Monte Carlo Uncertainty")
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


###############################################################################
# Plotting function 3
###############################################################################
def plot_idl_style(dem):
    """
    Faithful mirror of IDL's xrt_dem_iterative2.pro plotting style.

    - Black dotted lines -> Monte Carlo DEMs (if available)
    - Green line -> Best-fit DEM
    """
    if not hasattr(dem, "dem"):
        dem.solve()

    logT = dem.logT
    best_fit = dem.dem

    fig, ax = plt.subplots(figsize=(8, 6))

    if hasattr(dem, "_dem_ensemble"):
        mc_dems = np.array(dem._dem_ensemble)
        for i in range(mc_dems.shape[0]):
            ax.step(
                logT,
                np.log10(mc_dems[i] + 1e-40),
                where="mid",
                linestyle=":",
                color="black",
                alpha=0.3,
                linewidth=0.6,
            )

    ax.step(
        logT,
        np.log10(best_fit + 1e-40),
        where="mid",
        color="green",
        linewidth=2,
        label="Best-fit DEM",
    )

    ax.set_xlabel("log10 T [K]")
    ax.set_ylabel("log10 DEM [cm$^{-5}$ K$^{-1}$]")
    ax.set_xlim(logT.min(), logT.max())
    ax.set_ylim(
        np.floor(np.log10(best_fit.min() + 1e-40)),
        np.ceil(np.log10(best_fit.max() + 1e-40)),
    )
    ax.set_title("DEM (IDL Style)")
    ax.legend()
    ax.grid(alpha=0.3)

    plt.tight_layout()
    plt.show()


###############################################################################
# Plotting function 4
###############################################################################
def plot_fit_residuals(dem):
    """
    Plot observed vs fitted intensities and residuals.
    """
    if not hasattr(dem, "dem"):
        dem.solve()

    if not hasattr(dem, "fitted_intensities"):
        raise AttributeError(
            "No fitted intensities found. Run fit_dem() or solve() first."
        )

    obs = dem._observed_intensities
    fit = dem.fitted_intensities
    sigma = dem.intensity_errors.to_value(u.DN / u.s)

    filters = dem.filter_names
    indices = np.arange(len(obs))

    plt.figure(figsize=(7, 5))
    plt.errorbar(indices, obs, yerr=sigma, fmt="o", label="Observed", color="black")
    plt.plot(indices, fit, "s", label="Fitted", color="red")
    plt.xticks(indices, filters, rotation=45)
    plt.ylabel("Intensity [DN/s/pix]")
    plt.title("Observed vs Fitted Intensities")
    plt.legend()
    plt.tight_layout()
    plt.show()

    residuals = (obs - fit) / sigma
    plt.figure(figsize=(7, 4))
    plt.axhline(0, color="gray", linestyle="--")
    plt.plot(indices, residuals, "o", color="blue")
    plt.xticks(indices, filters, rotation=45)
    plt.ylabel("(Obs - Fit) / σ")
    plt.title("Residuals per Filter")
    plt.tight_layout()
    plt.show()


###############################################################################
# Plotting function 5
###############################################################################
def plot_dem_with_median_bins(dem):
    """
    Reproduce IDL-style DEM plot with:
    - Best-fit DEM (green)
    - Monte Carlo DEMs (dotted black)
    - Median DEM (blue)
    - Closest DEM to the median (orange)
    """
    if not hasattr(dem, "dem"):
        dem.solve()
    if not hasattr(dem, "_dem_ensemble"):
        raise AttributeError(
            "Monte Carlo ensemble not available. Run with monte_carlo_runs > 0."
        )

    logT = dem.logT
    mc_dems = np.array(dem._dem_ensemble)
    best_fit = dem.dem

    med = np.median(mc_dems, axis=0)
    diffs = np.linalg.norm(mc_dems - med, axis=1)
    closest_idx = np.argmin(diffs)
    closest_dem = mc_dems[closest_idx]

    fig, ax = plt.subplots(figsize=(9, 6))

    for i in range(mc_dems.shape[0]):
        ax.step(
            logT,
            np.log10(mc_dems[i] + 1e-40),
            where="mid",
            linestyle=":",
            color="black",
            alpha=0.3,
            linewidth=0.6,
        )

    ax.step(
        logT,
        np.log10(best_fit + 1e-40),
        where="mid",
        color="green",
        linewidth=2,
        label="Obs DEM",
    )
    ax.step(
        logT,
        np.log10(med + 1e-40),
        where="mid",
        color="blue",
        linewidth=1.8,
        label="Median in bins",
    )
    ax.step(
        logT,
        np.log10(closest_dem + 1e-40),
        where="mid",
        color="orange",
        linewidth=1.8,
        label="Closest DEM to median",
    )

    ax.set_xlim(dem.min_T, dem.max_T)
    ax.set_ylim(0, 30)
    ax.set_xlabel("Log T (K)")
    ax.set_ylabel("Log DEM [cm$^{-5}$ K$^{-1}$]")
    ax.legend()
    ax.grid(alpha=0.3)
    ax.set_title("DEM with Monte Carlo Spread, Median, and Closest Fit (IDL Style)")

    plt.tight_layout()
    plt.show()


###############################################################################
# Plotting function 6
###############################################################################
def plot_iteration_stats(dem):
    """
    Plot x^2 convergence across solver iterations.
    """
    if not hasattr(dem, "_iteration_chi2") or len(dem._iteration_chi2) == 0:
        raise AttributeError(
            "No iteration stats found. Run fit_dem() or solve() first."
        )

    chi2_vals = np.array(dem._iteration_chi2)

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(range(len(chi2_vals)), chi2_vals, lw=1.5)
    ax.set_xlabel("Iteration")
    ax.set_ylabel("Chi²")
    ax.set_title("Chi² Convergence")
    ax.grid(alpha=0.3)

    if chi2_vals.max() / max(chi2_vals.min(), 1e-10) > 1e4:
        ax.set_yscale("log")

    plt.tight_layout()
    plt.show()

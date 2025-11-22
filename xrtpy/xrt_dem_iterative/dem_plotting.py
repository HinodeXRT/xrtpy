__all__ = [
    "plot_dem",
    "plot_dem_mc",
    "plot_observed_vs_modeled",
]

import astropy.units as u
import matplotlib.pyplot as plt
import numpy as np


def plot_dem(solver):
    """
    Plot the base DEM solution (no Monte Carlo needed).

    Parameters
    ----------
    solver : XRTDEMIterative
        A solved DEM object (xrtpy.xrt_dem_iterative.XRTDEMIterative).
        If .dem does not exist yet, solve() will be called.
    """
    # Ensure base DEM is available
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    dem = np.asarray(solver.dem, dtype=float)

    # Avoid log10(0)
    dem_safe = np.clip(dem, 1e-40, None)

    plt.figure(figsize=(8, 5))
    plt.step(
        logT,
        np.log10(dem_safe),
        where="mid",
        color="black",
        linewidth=2.0,
        label="Base DEM",
    )

    plt.xlabel(r"log$_{10} T$  [K]")
    plt.ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    plt.title("Base DEM Solution")
    plt.grid(True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # Nice y-limits based on finite values
    finite = np.isfinite(np.log10(dem_safe))
    if np.any(finite):
        y = np.log10(dem_safe[finite])
        pad = 0.3 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_dem_mc(
    solver,
    mc_color="gray",
    base_color="red",
    alpha_mc=0.2,
    lw_mc=1.0,
    lw_base=2.0,
    show_envelope=True,
):
    """
    Plot DEM with Monte Carlo ensemble.

    - Base DEM in red.
    - Each Monte Carlo realization as a thin gray step curve.
    - Optional 68% envelope (16–84 percentile) if MC is available.

    If no Monte Carlo results are present, this gracefully falls back
    to plotting only the base DEM.
    """
    # Ensure base DEM is available
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    base_dem = np.asarray(solver.dem, dtype=float)
    base_dem_safe = np.clip(base_dem, 1e-40, None)

    has_mc = hasattr(solver, "mc_dem") and solver.mc_dem is not None

    if has_mc:
        mc_dem = np.asarray(solver.mc_dem, dtype=float)  # shape (N+1, n_T)
        # If someone set monte_carlo_runs=0, mc_dem may be (1, n_T)
        N = max(0, mc_dem.shape[0] - 1)
    else:
        mc_dem = None
        N = 0

    plt.figure(figsize=(9, 6))

    # --- Monte Carlo curves (index 1..N) ---
    if has_mc and N > 0:
        for i in range(1, N + 1):
            dem_i = np.clip(mc_dem[i, :], 1e-40, None)
            plt.step(
                logT,
                np.log10(dem_i),
                where="mid",
                color=mc_color,
                alpha=alpha_mc,
                linewidth=lw_mc,
            )

        # Optional envelope (16–84 percentile over MC curves only)
        if show_envelope and N > 1:
            dem_low = np.percentile(mc_dem[1:, :], 16, axis=0)
            dem_high = np.percentile(mc_dem[1:, :], 84, axis=0)

            dem_low = np.clip(dem_low, 1e-40, None)
            dem_high = np.clip(dem_high, 1e-40, None)

            plt.fill_between(
                logT,
                np.log10(dem_low),
                np.log10(dem_high),
                step="mid",
                color=mc_color,
                alpha=0.2,
                label="68% interval (MC)",
            )

    # --- Base DEM ---
    plt.step(
        logT,
        np.log10(base_dem_safe),
        where="mid",
        color=base_color,
        linewidth=lw_base,
        label="Base DEM",
    )

    plt.xlabel(r"log$_{10} T$  [K]")
    plt.ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    title_suffix = f" ({N} MC realizations)" if has_mc and N > 0 else ""
    plt.title("DEM with Monte Carlo" + title_suffix)
    plt.grid(True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # y-limits based on base DEM and MC if available
    logs = [np.log10(base_dem_safe)]
    if has_mc and N > 0:
        logs.append(np.log10(np.clip(mc_dem[1:, :], 1e-40, None)).ravel())
    logs_all = np.concatenate(logs)
    finite = np.isfinite(logs_all)
    if np.any(finite):
        y = logs_all[finite]
        pad = 0.3 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.legend()
    plt.tight_layout()
    plt.show()

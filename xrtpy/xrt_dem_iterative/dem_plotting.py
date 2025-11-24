__all__ = [
    "plot_dem",
    "plot_dem_mc",
]

import matplotlib.pyplot as plt
import numpy as np


def plot_dem(solver):
    """
    Plot the base DEM solution in log10 psace (no Monte Carlo needed).

    Parameters
    ----------
    solver : XRTDEMIterative
        Fully initialized DEM solver.
    """
    # Ensure base DEM is available
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    dem = np.asarray(solver.dem, dtype=float)

    # Avoid log10(0)
    dem_safe = np.clip(dem, 1e-40, None)
    log_dem = np.log10(dem_safe)

    plt.figure(figsize=(8, 5))
    plt.step(
        logT,
        log_dem,
        where="mid",
        color="black",
        linewidth=2.0,
        label="Base DEM",
    )

    plt.xlabel(r"log$_{10} T$  [K]")
    plt.ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    plt.title("Base DEM Solution")
    plt.grid(visible=True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # y-limits based on finite values
    finite = np.isfinite(log_dem)
    if np.any(finite):
        y = np.log10(dem_safe[finite])
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_dem_mc(
    solver,
    mc_color="black",
    base_color="#1E90FF",
    alpha_mc=0.18,  # 0.2
    lw_mc=1.0,
    lw_base=2.0,
    figsize=(9, 6),
):
    """
    Plot DEM with Monte Carlo ensemble (if present).

    - Base DEM: thick colored step curve - Dodger Blue
    - MC curves: thin transparent Black
    - Automatically chooses limits, even if MC not present

    If no Monte Carlo results are present, this gracefully falls back
    to plotting only the base DEM.
    """

    # Ensure base DEM is available
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    base_dem = np.asarray(solver.dem, dtype=float)
    base_safe_dem = np.clip(base_dem, 1e-100, None)  # 1e-40
    log_base_dem = np.log10(base_safe_dem)

    # CHecking for Monte Carlo
    has_mc = hasattr(solver, "mc_dem") and solver.mc_dem is not None
    if has_mc:
        mc_dem = np.asarray(solver.mc_dem, dtype=float)  # shape (N+1, n_T)
        # If someone set monte_carlo_runs=0, mc_dem may be (1, n_T)
        N = max(0, mc_dem.shape[0] - 1)
    else:
        mc_dem = None
        N = 0

    plt.figure(figsize=figsize)

    # Plot Monte Carlo curves (index 1..N)
    if has_mc and N > 0:
        for i in range(1, N + 1):
            dem_i = np.clip(mc_dem[i, :], 1e-100, None)  # 1e-40
            plt.step(
                logT,
                np.log10(dem_i),
                where="mid",
                color=mc_color,
                alpha=alpha_mc,
                linewidth=lw_mc,
            )

    # Base DEM
    plt.step(
        logT,
        log_base_dem,
        where="mid",
        color=base_color,
        linewidth=lw_base,
        label="Base DEM",
    )

    suffix = f" ({N} MC realizations)" if N > 0 else ""
    plt.title("DEM with Monte Carlo" + suffix)
    plt.xlabel(r"log$_{10} T$  [K]")
    plt.ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    plt.grid(visible=True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # y-limits based on base DEM and MC if available
    logs = [np.log10(log_base_dem)]

    if has_mc and N > 0:
        logs.append(np.log10(np.clip(mc_dem[1:, :], 1e-100, None)).ravel())  # 1e-40

    logs_all = np.concatenate(logs)

    finite = np.isfinite(logs_all)
    if np.any(finite):
        y = logs_all[finite]
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.legend()
    plt.tight_layout()
    plt.show()

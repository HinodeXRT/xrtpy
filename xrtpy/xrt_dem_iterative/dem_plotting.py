__all__ = [
    "plot_dem",
    "plot_dem_mc",
]

import matplotlib.pyplot as plt
import numpy as np


def plot_dem(solver):
    """
    Plot the base DEM solution in log10 space (no Monte Carlo needed).

    Parameters
    ----------
    solver : XRTDEMIterative
        Fully initialized DEM solver.
    """
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    dem = np.asarray(solver.dem, dtype=float)

    dem_safe = np.clip(dem, 1e-40, None)
    log_dem = np.log10(dem_safe)

    # Build title from filter names if available
    filters = getattr(solver, "filter_names", None)
    if filters:
        filter_str = ", ".join(filters)
        title = f"Hinode/XRT — Base DEM Solution\nFilters: {filter_str}"
    else:
        title = "Hinode/XRT — Base DEM Solution"

    plt.figure(figsize=(8, 5))
    plt.step(
        logT,
        log_dem,
        where="mid",
        color="black",
        linewidth=2.0,
        label="Base DEM",
    )

    plt.xlabel(r"$\log_{10}\,T$  [K]")
    plt.ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    plt.title(title)
    plt.grid(visible=True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # y-limits based on finite values
    finite = np.isfinite(log_dem)
    if np.any(finite):
        y = log_dem[finite]
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.legend()
    plt.tight_layout()
    plt.show()


def plot_dem_mc(
    solver,
    mc_color="black",
    base_color="#1E90FF",
    alpha_mc=0.18,
    lw_mc=1.0,
    lw_base=2.0,
    figsize=(9, 6),
):
    """
    Plot DEM with Monte Carlo ensemble (if present).

    - Base DEM: thick Dodger Blue step curve
    - MC realizations: thin transparent black step curves
    - Falls back gracefully to base-only if no MC results are present
    """
    if not hasattr(solver, "dem"):
        solver.solve()

    logT = solver.logT
    base_dem = np.asarray(solver.dem, dtype=float)
    base_safe_dem = np.clip(base_dem, 1e-100, None)
    log_base_dem = np.log10(base_safe_dem)

    has_mc = hasattr(solver, "mc_dem") and solver.mc_dem is not None
    if has_mc:
        mc_dem = np.asarray(solver.mc_dem, dtype=float)  # shape (N+1, n_T)
        N = max(0, mc_dem.shape[0] - 1)
    else:
        mc_dem = None
        N = 0

    # Build title from filter names if available
    filters = getattr(solver, "filter_names", None)
    if filters:
        filter_str = ", ".join(filters)
        title = f"Hinode/XRT — DEM with Monte Carlo\nFilters: {filter_str}"
    else:
        title = "Hinode/XRT — DEM with Monte Carlo"

    if N > 0:
        title += f"  ({N} MC realizations)"

    plt.figure(figsize=figsize)

    # Plot MC curves first (so base DEM renders on top)
    mc_handle = None
    if has_mc and N > 0:
        for i in range(1, N + 1):
            dem_i = np.clip(mc_dem[i, :], 1e-100, None)
            (line,) = plt.step(
                logT,
                np.log10(dem_i),
                where="mid",
                color=mc_color,
                alpha=alpha_mc,
                linewidth=lw_mc,
            )
            if i == 1:
                mc_handle = line  # keep only one handle for the legend

    # Base DEM on top
    (base_handle,) = plt.step(
        logT,
        log_base_dem,
        where="mid",
        color=base_color,
        linewidth=lw_base,
        label="Base DEM",
    )

    # Legend: always show base DEM; add MC entry only if MC curves are present
    if has_mc and N > 0 and mc_handle is not None:
        mc_handle.set_alpha(1.0)  # full opacity for the legend swatch
        plt.legend(
            handles=[base_handle, mc_handle],
            labels=["Base DEM", f"MC realizations (N={N})"],
        )
    else:
        plt.legend(handles=[base_handle], labels=["Base DEM"])

    plt.title(title)
    plt.xlabel(r"$\log_{10}\,T$  [K]")
    plt.ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    plt.grid(visible=True, alpha=0.3)
    plt.xlim(logT.min(), logT.max())

    # y-limits — bug fix: log_base_dem is already in log space, do not take log10 again
    logs = [log_base_dem]
    if has_mc and N > 0:
        logs.append(np.log10(np.clip(mc_dem[1:, :], 1e-100, None)).ravel())
    logs_all = np.concatenate(logs)

    finite = np.isfinite(logs_all)
    if np.any(finite):
        y = logs_all[finite]
        pad = 0.25 * (y.max() - y.min() + 1e-6)
        plt.ylim(y.min() - pad, y.max() + pad)

    plt.tight_layout()
    plt.show()

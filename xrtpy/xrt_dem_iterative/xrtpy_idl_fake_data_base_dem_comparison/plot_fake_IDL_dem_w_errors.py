from pathlib import Path

import matplotlib
import numpy as np

# No GUI needed since we're only saving PNGs
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from scipy.io import readsav

# -----------------------------------------------------------------------------
# Directories
# -----------------------------------------------------------------------------

BASE_DIR = Path(
    "/Users/jvelasq/Projects/xrtpy/xrtpy/xrt_dem_iterative/"
    "xrtpy_idl_fake_data_base_dem_comparison"
)

IDL_DIR = BASE_DIR / "IDL_fake_DEM_data_w_errors"
OUT_DIR = BASE_DIR / "IDL_fake_DEM_plots"

OUT_DIR.mkdir(parents=True, exist_ok=True)


# -----------------------------------------------------------------------------
# Load IDL .sav file
# -----------------------------------------------------------------------------


def load_idl_base(dem_id, nfilters):
    """
    Load IDL base DEM .sav file.

    Returns
    -------
    logT : ndarray, shape (26,)
    dem : ndarray, shape (26,)
    chisq : float
    """
    sav_file = (
        IDL_DIR / f"idl_dem_idx{dem_id}_base_wErrors_10Oct2008_{nfilters}filters.sav"
    )

    if not sav_file.exists():
        raise FileNotFoundError(f"Cannot find:\n{sav_file}")

    data = readsav(str(sav_file), python_dict=True)

    logT = np.asarray(data["logt_out"]).ravel()
    dem = np.asarray(data["dem_out"]).ravel()
    chisq = float(np.asarray(data["chisq"]).ravel()[0])

    return logT, dem, chisq


# -----------------------------------------------------------------------------
# Plot one DEM
# -----------------------------------------------------------------------------


def plot_idl_dem(dem_id, nfilters):

    logT, dem, chisq = load_idl_base(dem_id, nfilters)

    # Avoid log10(0)
    log_dem = np.log10(np.maximum(dem, 1e-40))

    fig, ax = plt.subplots(figsize=(8, 5))

    ax.step(
        logT,
        log_dem,
        where="mid",
        color="orange",
        linewidth=2.5,
        label=rf"IDL  |  $\chi^2$ = {chisq:.4f}",
    )

    ax.set_title(f"IDL Base DEM\nDEM Index {dem_id} — {nfilters} Filters")

    ax.set_xlabel(r"$\log_{10}(T)$ [K]")
    ax.set_ylabel(r"$\log_{10}(\mathrm{DEM})$ [cm$^{-5}$ K$^{-1}$]")

    ax.set_xlim(5.5, 8.0)
    ax.set_ylim(16, 26)

    ax.grid(True, alpha=0.3)
    ax.legend()

    fig.tight_layout()

    outfile = (
        OUT_DIR / f"idl_dem_idx{dem_id}_base_wErrors_10Oct2008_{nfilters}filters.png"
    )

    fig.savefig(
        outfile,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)

    print(f"Saved: {outfile.name}")


# -----------------------------------------------------------------------------
# Make all 24 plots
# -----------------------------------------------------------------------------

if __name__ == "__main__":
    dem_ids = [14, 15, 65, 99, 112, 159]
    filter_counts = [4, 5, 6, 7]

    for dem_id in dem_ids:
        for nfilters in filter_counts:
            try:
                plot_idl_dem(dem_id, nfilters)
            except Exception as e:
                print(f"FAILED: DEM {dem_id}, {nfilters} filters\n{e}\n")

    print(f"\nAll plots saved to:\n{OUT_DIR}")

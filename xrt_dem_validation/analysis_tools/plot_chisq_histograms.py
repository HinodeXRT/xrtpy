"""
Chi-square histograms for the 162-DEM synthetic validation campaign.

Two overlaid figures (IDL orange vs XRTpy blue), log10(chi2) x-axis:

  1. BASE      -- per-DEM base chi-square (one value per DEM per solver)
                  plots/histograms/chisq_hist_base.png
  2. MC-MEDIAN -- per-DEM median of the 1000 MC chi-square values
                  plots/histograms/chisq_hist_mc_median.png

Each figure marks the stall threshold (chi2 = 100) with a vertical line
and reports, per solver, how many DEMs fall beyond it.

Data read directly from source:
  IDL base : data/all_IDL_synthetic_dem_data/*_base.sav          (chisq[0])
  IDL MC   : data/all_IDL_MC_synthetic_dem_data/*_MC1000.sav     (chisq array)
  XRTpy    : xrtpy_output/*_MC1000.npz                           (mc_chisq)

Pure load-and-plot. Uses only DEMs present for BOTH solvers.

Usage (from synthetic_dem_validation/):
    python plot_chisq_histograms.py
"""

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
IDL_MC_DIR = BASE_DIR / "data" / "all_IDL_MC_synthetic_dem_data"
XRTPY_DIR = BASE_DIR / "xrtpy_output"
PLOT_DIR = BASE_DIR / "plots" / "histograms"

N_MC = 1000
STALL_THRESHOLD = 100.0
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"


# ── Helpers ───────────────────────────────────────────────────────────────────
def _dem_id(p: Path) -> int:
    m = re.search(r"idx(\d+)_", p.name)
    return int(m.group(1)) if m else -1


def idl_base_chisq(dem_id):
    path = IDL_BASE_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"
    return float(np.atleast_1d(readsav(str(path), python_dict=True)["chisq"])[0])


def idl_mc_chisq(dem_id):
    """Full MC chi-square array from the IDL MC .sav, dropping column 0."""
    path = IDL_MC_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_MC{N_MC}.sav"
    chi = np.atleast_1d(readsav(str(path), python_dict=True)["chisq"]).astype(float)
    return chi[1:] if chi.size > 1 else chi


def xrtpy_chisq(dem_id):
    """Returns (base_chisq, mc_chisq_array) from the XRTpy npz."""
    path = XRTPY_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_MC{N_MC}.npz"
    d = np.load(path, allow_pickle=True)
    chi = np.asarray(d["mc_chisq"], dtype=float)
    return float(chi[0]), chi[1:]


# ── Histogram figure ──────────────────────────────────────────────────────────
def make_hist(idl_vals, py_vals, title, out_png, xlabel):
    idl_vals = np.asarray(idl_vals, dtype=float)
    py_vals = np.asarray(py_vals, dtype=float)

    # log10 axis; guard against non-positive
    log_idl = np.log10(np.clip(idl_vals, 1e-6, None))
    log_py = np.log10(np.clip(py_vals, 1e-6, None))

    lo = np.floor(min(log_idl.min(), log_py.min()))
    hi = np.ceil(max(log_idl.max(), log_py.max()))
    bins = np.linspace(lo, hi, 40)

    median_idl = float(np.median(idl_vals))
    median_py = float(np.median(py_vals))

    fig, ax = plt.subplots(figsize=(10, 6))

    ax.hist(log_idl, bins=bins, color=COLOR_IDL, alpha=0.55, label="IDL")
    ax.hist(log_py, bins=bins, color=COLOR_XRTPY, alpha=0.55, label="XRTpy")

    # Median chi-square drawn at log10(median); robust to the stall outliers
    ax.axvline(np.log10(max(median_idl, 1e-6)), color=COLOR_IDL,
               linestyle="--", linewidth=2.0,
               label=f"IDL median χ² = {median_idl:.1f}")
    ax.axvline(np.log10(max(median_py, 1e-6)), color=COLOR_XRTPY,
               linestyle="--", linewidth=2.0,
               label=f"XRTpy median χ² = {median_py:.1f}")

    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.set_ylabel("Number of DEMs", fontsize=12)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=10, loc="upper right")

    fig.tight_layout()
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    print(f"Saved: {out_png}")


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    idl_ids = {
        _dem_id(p)
        for p in IDL_MC_DIR.glob(f"idl_synthetic_dem_idx*_allfilters_MC{N_MC}.sav")
    }
    py_ids = {
        _dem_id(p)
        for p in XRTPY_DIR.glob(f"xrtpy_dem_idx*_allfilters_MC{N_MC}.npz")
    }
    dem_ids = sorted(idl_ids & py_ids)
    if not dem_ids:
        raise FileNotFoundError("No DEMs with both IDL MC and XRTpy MC present.")

    missing = sorted((idl_ids | py_ids) - (idl_ids & py_ids))
    if missing:
        print(f"Note: {len(missing)} DEM(s) skipped (one solver only): "
              f"{missing[:15]}{' ...' if len(missing) > 15 else ''}")
    print(f"Using {len(dem_ids)} DEMs present for both solvers.")

    idl_base, py_base = [], []
    idl_mcmed, py_mcmed = [], []

    for dem_id in dem_ids:
        idl_base.append(idl_base_chisq(dem_id))
        idl_mcmed.append(np.median(idl_mc_chisq(dem_id)))

        py_b, py_mc = xrtpy_chisq(dem_id)
        py_base.append(py_b)
        py_mcmed.append(np.median(py_mc))

    # Figure 1: base
    make_hist(
        idl_base, py_base,
        title=f"Base DEM χ² Distribution — {len(dem_ids)} synthetic DEMs",
        out_png=PLOT_DIR / "chisq_hist_base.png",
        xlabel=r"log$_{10}$($\chi^2$)  (base DEM)",
    )

    # Figure 2: MC median
    make_hist(
        idl_mcmed, py_mcmed,
        title=f"Median MC χ² Distribution — {len(dem_ids)} synthetic DEMs "
              f"({N_MC} MC each)",
        out_png=PLOT_DIR / "chisq_hist_mc_median.png",
        xlabel=r"log$_{10}$(median $\chi^2$)  (1000 MC per DEM)",
    )

    print("\nDone.")


if __name__ == "__main__":
    main()

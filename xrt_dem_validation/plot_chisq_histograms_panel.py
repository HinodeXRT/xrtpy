"""
Side-by-side 1x2 panel of the two chi-square histograms, for the poster.

Left  : Base DEM chi-square (162 IDL + 162 XRTpy)
Right : All-MC chi-square (~162,000 IDL + ~162,000 XRTpy)

Exactly the same content/style as plot_chisq_histograms_v2.py, just placed
in one figure. Each panel keeps its own axes (counts differ: 162 vs 162,000).
One shared legend across the top; per-panel median lines with values.

Sources (each solver from its own clean file):
  Python : data/xrtpy_monte_carlo_synthetic_dem_data/*_MC1000.npz  (mc_chisq[0] base, mc_chisq[1:] MC)
  IDL    : data/idl_base_chisq.csv, data/idl_mc_chisq_long.csv (mc_run>=1)

Usage (from synthetic_dem_validation/):
    python plot_chisq_histograms_panel_v2.py
"""

import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.patches import Patch

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
XRTPY_DIR = BASE_DIR / "data" / "xrtpy_monte_carlo_synthetic_dem_data"
IDL_BASE_CSV = BASE_DIR / "data" / "idl_base_chisq.csv"
IDL_MC_CSV = BASE_DIR / "data" / "idl_mc_chisq_long.csv"
PLOT_DIR = BASE_DIR / "plots" / "histograms"

N_MC = 1000
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"
CHI_FLOOR = 1e-6


# ── Loaders ───────────────────────────────────────────────────────────────────
def _dem_id(p):
    return int(re.search(r"idx(\d+)_", p.name).group(1))


def xrtpy_base_and_mc():
    base, mc = [], []
    files = sorted(XRTPY_DIR.glob(f"xrtpy_dem_idx*_allfilters_MC{N_MC}.npz"),
                   key=_dem_id)
    if not files:
        raise FileNotFoundError(f"No XRTpy MC npz in {XRTPY_DIR}")
    for p in files:
        d = np.load(p, allow_pickle=True)
        chi = np.asarray(d["mc_chisq"], dtype=float)
        base.append(chi[0])
        mc.append(chi[1:])
    return np.array(base), np.concatenate(mc)


def idl_base():
    vals = []
    with open(IDL_BASE_CSV) as f:
        for row in csv.DictReader(f):
            vals.append(float(row["chisq"]))
    return np.array(vals)


def idl_mc():
    vals = []
    with open(IDL_MC_CSV) as f:
        reader = csv.reader(f)
        next(reader)
        for dem_id, mc_run, chisq in reader:
            if int(mc_run) >= 1:
                vals.append(float(chisq))
    return np.array(vals)


# ── One panel ─────────────────────────────────────────────────────────────────
def draw_panel(ax, idl_vals, py_vals, title, xlabel):
    idl_vals = np.asarray(idl_vals, dtype=float)
    py_vals = np.asarray(py_vals, dtype=float)

    log_idl = np.log10(np.clip(idl_vals, CHI_FLOOR, None))
    log_py = np.log10(np.clip(py_vals, CHI_FLOOR, None))

    lo = np.floor(min(log_idl.min(), log_py.min()))
    hi = np.ceil(max(log_idl.max(), log_py.max()))
    bins = np.linspace(lo, hi, 60)

    med_idl = float(np.median(idl_vals))
    med_py = float(np.median(py_vals))

    ax.hist(log_idl, bins=bins, color=COLOR_IDL, alpha=0.55)
    ax.hist(log_py, bins=bins, color=COLOR_XRTPY, alpha=0.55)
    ax.axvline(np.log10(max(med_idl, CHI_FLOOR)), color=COLOR_IDL,
               linestyle="--", linewidth=2.0)
    ax.axvline(np.log10(max(med_py, CHI_FLOOR)), color=COLOR_XRTPY,
               linestyle="--", linewidth=2.0)

    # Per-panel median annotation (since medians differ between panels)
    ax.text(0.97, 0.80,
            f"IDL median χ² = {med_idl:.2f}\nXRTpy median χ² = {med_py:.2f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=10,
            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))

    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel, fontsize=12)
    ax.grid(visible=True, alpha=0.3)


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    print("Loading XRTpy from npz, IDL from CSV...")
    py_base, py_mc = xrtpy_base_and_mc()
    idl_b = idl_base()
    idl_m = idl_mc()

    fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(18, 6.5))

    draw_panel(
        ax_left, idl_b, py_base,
        title=f"Base DEM χ²  —  {len(idl_b)} synthetic DEMs",
        xlabel=r"log$_{10}$($\chi^2$)  (base DEM)",
    )
    ax_left.set_ylabel("Count", fontsize=12)

    draw_panel(
        ax_right, idl_m, py_mc,
        title=f"All-MC χ²  —  162 DEMs × 1000 MC  "
              f"({len(idl_m):,} draws each)",
        xlabel=r"log$_{10}$($\chi^2$)  (every MC draw)",
    )
    ax_right.set_ylabel("Count", fontsize=12)

    # One shared legend across the top
    handles = [
        Patch(fc=COLOR_IDL, alpha=0.55, label="IDL"),
        Patch(fc=COLOR_XRTPY, alpha=0.55, label="XRTpy"),
        Line2D([], [], color=COLOR_IDL, ls="--", lw=2, label="IDL median"),
        Line2D([], [], color=COLOR_XRTPY, ls="--", lw=2, label="XRTpy median"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=12,
               frameon=True, bbox_to_anchor=(0.5, 1.0))

    fig.tight_layout(rect=[0, 0, 1, 0.94])

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = PLOT_DIR / "chisq_hist_panel_base_and_mc_v2.png"
    fig.savefig(out_png, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()

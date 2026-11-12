"""
Poster 2x2 grid — BASE DEMs (synthetic + IDL + XRTpy), no Monte Carlo.

Give four DEM indices on the command line; they fill the grid in the order
you pass them (top-left, top-right, bottom-left, bottom-right). Shared axes:
one x-label along the bottom, one y-label up the side, per-panel legends
kept compact. Contribution bands shaded if data/contribution_bands.csv exists.

Loads per DEM:
    - Band:    data/contribution_bands.csv (optional)
    - True:    data/synthetic_dems_data/DEM_{id}.txt
    - IDL:     data/IDL_synthetic_dem_data/*_base.sav
    - XRTpy:   data/xrtpy_monte_carlo_synthetic_dem_data/*_MC1000.npz (row 0) or *_base.npz fallback

Saved to plots/poster/poster_base_{a}_{b}_{c}_{d}.png

Usage (from synthetic_dem_validation/):
    python poster_grid_base_v2.py 14 65 112 159
    python poster_grid_base_v2.py 14 65 112 159 --out my_panel.png
"""

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
IDL_BASE_DIR = BASE_DIR / "data" / "IDL_synthetic_dem_data"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
XRTPY_DIR = BASE_DIR / "data" / "xrtpy_monte_carlo_synthetic_dem_data"
BANDS_CSV = BASE_DIR / "data" / "contribution_bands.csv"
PLOT_DIR = BASE_DIR / "plots" / "poster"

COLOR_TRUE = "red"
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"
BAND_COLOR = "0.6"
BAND_ALPHA = 0.25
GRID_STEP = 0.1
YLIM = (14, 26)
XLIM = (5.5, 8.0)


# ── Helpers ───────────────────────────────────────────────────────────────────
def _log10_dem(dem, floor=1e-40):
    return np.log10(np.maximum(np.asarray(dem, dtype=float), floor))


def _runs_by_T(dem, nT):
    dem = np.atleast_2d(np.array(dem))
    if dem.shape[1] == nT:
        return dem
    if dem.shape[0] == nT:
        return dem.T
    raise ValueError(f"Cannot orient DEM array {dem.shape} with nT={nT}")


def _decode(arr):
    return [x.decode() if isinstance(x, bytes) else str(x)
            for x in np.atleast_1d(arr)]


def load_bands():
    bands = {}
    if not BANDS_CSV.exists():
        return bands
    with open(BANDS_CSV) as f:
        for row in csv.DictReader(f):
            temps = [float(x) for x in row["temps"].split()] if row["temps"] else []
            bands[int(row["dem_id"])] = temps
    return bands


def contiguous_spans(temps, step=GRID_STEP):
    if not temps:
        return []
    t = np.sort(np.asarray(temps, dtype=float))
    spans, start, prev = [], t[0], t[0]
    for x in t[1:]:
        if x - prev > step * 1.5:
            spans.append((start, prev)); start = x
        prev = x
    spans.append((start, prev))
    half = step / 2.0
    return [(lo - half, hi + half) for lo, hi in spans]


def load_idl_base(dem_id):
    path = IDL_BASE_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"
    data = readsav(str(path), python_dict=True)
    logT = np.array(data["logt_out"]).ravel()
    dem = _runs_by_T(data["dem_out"], len(logT))[0]
    chisq = float(np.atleast_1d(data["chisq"])[0])
    return logT, dem, chisq


def load_xrtpy_base(dem_id):
    mc = XRTPY_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_MC1000.npz"
    base = XRTPY_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
    if mc.exists():
        d = np.load(mc, allow_pickle=True)
        return d["logT"], d["mc_dem"][0], float(d["mc_chisq"][0])
    if base.exists():
        d = np.load(base, allow_pickle=True)
        return d["logT"], d["dem"], float(d["chisq"])
    raise FileNotFoundError(f"No XRTpy npz for DEM {dem_id}")


def load_true_dem(dem_id):
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    if not path.exists():
        return None
    d = np.loadtxt(path)
    return d[:, 0], _log10_dem(d[:, 1])


# ── Panel ─────────────────────────────────────────────────────────────────────
def draw_panel(ax, dem_id, bands):
    logT, idl_base, idl_chisq = load_idl_base(dem_id)
    _, py_base, py_chisq = load_xrtpy_base(dem_id)

    for lo, hi in contiguous_spans(bands.get(dem_id, [])):
        ax.axvspan(lo, hi, color=BAND_COLOR, alpha=BAND_ALPHA)

    true_dem = load_true_dem(dem_id)
    if true_dem is not None:
        ax.step(true_dem[0], true_dem[1], where="mid",
                color=COLOR_TRUE, linewidth=2.2)
    ax.step(logT, _log10_dem(idl_base), where="mid",
            color=COLOR_IDL, linewidth=2.2)
    ax.step(logT, _log10_dem(py_base), where="mid",
            color=COLOR_XRTPY, linewidth=2.2, linestyle="--")

    ax.set_xlim(*XLIM)
    ax.set_ylim(*YLIM)
    ax.grid(visible=True, alpha=0.3)

    # Per-panel identifier + chi-square, compact, top-right
    txt = (f"DEM {dem_id}\n"
           f"IDL χ²={idl_chisq:.2f}\n"
           f"XRTpy χ²={py_chisq:.2f}")
    ax.text(0.97, 0.95, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", ec="0.7",
                                  alpha=0.85))


# ── Main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(description="Poster 2x2 base DEM grid")
    parser.add_argument("dems", type=int, nargs=4,
                        help="Four DEM indices (top-left, top-right, "
                             "bottom-left, bottom-right).")
    parser.add_argument("--out", type=str, default=None,
                        help="Output filename (into plots/poster/).")
    args = parser.parse_args()

    bands = load_bands()

    fig, axes = plt.subplots(2, 2, figsize=(11, 8.5), sharex=True, sharey=True)
    for ax, dem_id in zip(axes.flat, args.dems):
        draw_panel(ax, dem_id, bands)

    # One x-label bottom, one y-label side
    for ax in axes[-1, :]:
        ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=13)
    for ax in axes[:, 0]:
        ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=13)

    # Single shared legend across the top
    handles = [
        plt.Line2D([], [], color=COLOR_TRUE, lw=2.2, label="Synthetic DEM"),
        plt.Line2D([], [], color=COLOR_IDL, lw=2.2, label="IDL base"),
        plt.Line2D([], [], color=COLOR_XRTPY, lw=2.2, ls="--", label="XRTpy base"),
        plt.Rectangle((0, 0), 1, 1, fc=BAND_COLOR, alpha=BAND_ALPHA,
                      label=">2% contribution"),
    ]
    fig.legend(handles=handles, loc="upper center", ncol=4, fontsize=11,
               frameon=True, bbox_to_anchor=(0.5, 1.0))

    fig.tight_layout(rect=[0, 0, 1, 0.955])

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    if args.out:
        out_png = PLOT_DIR / args.out
    else:
        a, b, c, d = args.dems
        out_png = PLOT_DIR / f"poster_base_{a}_{b}_{c}_{d}.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)
    print(f"Saved: {out_png}")


if __name__ == "__main__":
    main()

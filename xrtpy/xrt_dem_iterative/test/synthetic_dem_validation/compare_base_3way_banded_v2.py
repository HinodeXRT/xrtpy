"""
Three-way BASE DEM comparison WITH contribution-band shading, for the
162-DEM synthetic validation campaign.

Synthetic truth (red) + IDL base (orange) + XRTpy base (blue dashed),
with grey shaded regions marking the temperatures that carry > N% of the
emission -- where the filters actually constrain the DEM. Each contiguous
run of qualifying logT bins (spacing ~0.1) is shaded as its own box, so
multi-component DEMs get multiple bands with gaps between them.

Loads:
    - Bands:     data/contribution_bands.csv (from extract_contribution_bands.py)
    - True DEM:  data/synthetic_dems_data/DEM_{id}.txt
    - IDL base:  data/all_IDL_synthetic_dem_data/*_base.sav
    - XRTpy:     xrtpy_output/*_MC1000.npz (row 0), or *_base.npz fallback

Saved to plots/compare_base_banded/compare_base_dem_idx{id}_allfilters.png

Usage (from synthetic_dem_validation/):
    python compare_synthetic_base_3way_banded.py --dem 137   # one DEM
    python compare_synthetic_base_3way_banded.py             # all found
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
XRTPY_DIR = BASE_DIR / "xrtpy_output"
BANDS_CSV = BASE_DIR / "data" / "contribution_bands.csv"
PLOT_DIR = BASE_DIR / "plots" / "compare_base_banded"

COLOR_TRUE = "red"
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"
BAND_COLOR = "0.6"
BAND_ALPHA = 0.25
GRID_STEP = 0.1


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
    out = []
    for x in np.atleast_1d(arr):
        out.append(x.decode() if isinstance(x, bytes) else str(x))
    return out


def load_bands():
    bands = {}
    if not BANDS_CSV.exists():
        print(f"Note: {BANDS_CSV} not found -- plotting without shading.")
        return bands
    with open(BANDS_CSV) as f:
        for row in csv.DictReader(f):
            temps = [float(x) for x in row["temps"].split()] if row["temps"] else []
            bands[int(row["dem_id"])] = temps
    return bands


def contiguous_spans(temps, step=GRID_STEP):
    """Group sorted logT values into contiguous [lo, hi] spans; gaps split."""
    if not temps:
        return []
    t = np.sort(np.asarray(temps, dtype=float))
    spans = []
    start = prev = t[0]
    for x in t[1:]:
        if x - prev > step * 1.5:
            spans.append((start, prev))
            start = x
        prev = x
    spans.append((start, prev))
    half = step / 2.0
    return [(lo - half, hi + half) for lo, hi in spans]


# ── Loaders ───────────────────────────────────────────────────────────────────
def load_idl_base(dem_id):
    path = IDL_BASE_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"
    data = readsav(str(path), python_dict=True)
    logT = np.array(data["logt_out"]).ravel()
    dem = _runs_by_T(data["dem_out"], len(logT))[0]
    chisq = float(np.atleast_1d(data["chisq"])[0])
    filters = _decode(data["obs_index"])
    obs_val = np.atleast_1d(data["obs_val"]).astype(float)
    obs_err = np.atleast_1d(data["obs_err"]).astype(float)
    return logT, dem, chisq, filters, obs_val, obs_err


def load_xrtpy_base(dem_id):
    mc = XRTPY_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_MC1000.npz"
    base = XRTPY_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
    if mc.exists():
        d = np.load(mc, allow_pickle=True)
        return d["logT"], d["mc_dem"][0], float(d["mc_chisq"][0])
    if base.exists():
        d = np.load(base, allow_pickle=True)
        return d["logT"], d["dem"], float(d["chisq"])
    raise FileNotFoundError(f"No XRTpy npz for DEM {dem_id} in {XRTPY_DIR}")


def load_true_dem(dem_id):
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    if not path.exists():
        return None
    d = np.loadtxt(path)
    return d[:, 0], _log10_dem(d[:, 1])


# ── Plot ──────────────────────────────────────────────────────────────────────
def plot_one_dem(dem_id, bands):
    idl_logT, idl_base, idl_chisq, filters, obs_val, obs_err = load_idl_base(dem_id)
    py_logT, py_base, py_chisq = load_xrtpy_base(dem_id)

    if not np.allclose(idl_logT, py_logT, atol=1e-6):
        raise ValueError(f"IDL and XRTpy logT grids differ for DEM {dem_id}!")
    logT = py_logT

    fig, ax = plt.subplots(figsize=(11, 7))

    spans = contiguous_spans(bands.get(dem_id, []))
    for k, (lo, hi) in enumerate(spans):
        ax.axvspan(lo, hi, color=BAND_COLOR, alpha=BAND_ALPHA,
                   label=">2% contribution" if k == 0 else None)

    true_dem = load_true_dem(dem_id)
    if true_dem is not None:
        ax.step(true_dem[0], true_dem[1], where="mid",
                color=COLOR_TRUE, linewidth=2.5, label="Synthetic DEM")

    ax.step(logT, _log10_dem(idl_base), where="mid",
            color=COLOR_IDL, linewidth=2.5,
            label=f"IDL base  |  χ² = {idl_chisq:.4f}")
    ax.step(logT, _log10_dem(py_base), where="mid",
            color=COLOR_XRTPY, linewidth=2.5, linestyle="--",
            label=f"XRTpy base  |  χ² = {py_chisq:.4f}")

    ax.set_title(
        f"Base DEM Comparison — IDL vs XRTpy vs Synthetic — DEM Index {dem_id} "
        f"(all {len(filters)} filters)",
        fontsize=13, pad=34,
    )

    entries = [
        f"{f} [{i:.4g} ± {e:.3g}]"
        for f, i, e in zip(filters, obs_val, obs_err, strict=True)
    ]
    per_line = 4
    filter_label = "\n".join(
        ",   ".join(entries[k : k + per_line])
        for k in range(0, len(entries), per_line)
    )
    fig.text(0.5, 0.945, filter_label, ha="center", va="top", fontsize=7.5)

    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)
    ax.set_xlim(logT.min(), logT.max())
    ax.set_ylim(18, 26)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.87])
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = PLOT_DIR / f"compare_base_dem_idx{dem_id}_allfilters.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return idl_chisq, py_chisq, len(spans)


# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(p: Path) -> int:
    m = re.search(r"idx(\d+)_", p.name)
    return int(m.group(1)) if m else -1


def main():
    parser = argparse.ArgumentParser(description="3-way base comparison (banded)")
    parser.add_argument("--dem", type=int, default=None,
                        help="Single DEM index (e.g., 137). Omit for all found.")
    args = parser.parse_args()

    bands = load_bands()

    sav_files = sorted(
        IDL_BASE_DIR.glob("idl_synthetic_dem_idx*_allfilters_base.sav"),
        key=_dem_id_from_path,
    )
    if not sav_files:
        raise FileNotFoundError(f"No IDL base .sav files found in {IDL_BASE_DIR}")

    dem_ids = [_dem_id_from_path(p) for p in sav_files]
    if args.dem is not None:
        if args.dem not in dem_ids:
            raise FileNotFoundError(f"No IDL base .sav for DEM {args.dem}")
        dem_ids = [args.dem]

    n_no_band = 0
    print(f"Plotting {len(dem_ids)} DEM(s).")
    for k, dem_id in enumerate(dem_ids, 1):
        idl_chi, py_chi, n_spans = plot_one_dem(dem_id, bands)
        if n_spans == 0:
            n_no_band += 1
        print(f"[{k:3d}/{len(dem_ids)}] DEM {dem_id:3d}  "
              f"IDL {idl_chi:12.4f}   XRTpy {py_chi:10.4f}   bands: {n_spans}")

    if n_no_band:
        print(f"\nNote: {n_no_band} DEM(s) plotted without shading "
              f"(not in {BANDS_CSV.name}).")
    print("\nDone.")


if __name__ == "__main__":
    main()
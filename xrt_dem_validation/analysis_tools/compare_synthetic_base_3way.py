"""
Three-way BASE DEM comparison for the 162-DEM synthetic validation campaign:
synthetic truth + IDL base + XRTpy base (no Monte Carlo).

For each DEM index, loads:
    - True DEM:  data/synthetic_dems_data/DEM_{id}.txt
    - IDL base:  data/all_IDL_synthetic_dem_data/
                     idl_synthetic_dem_idx{id}_allfilters_base.sav
    - XRTpy:     xrtpy_output/xrtpy_dem_idx{id}_allfilters_MC1000.npz
                 (falls back to ..._base.npz if the MC file is absent)

and produces one plot per DEM:
    - True DEM   : red solid
    - IDL base   : orange solid   (chi-square in legend)
    - XRTpy base : dodger blue dashed (chi-square in legend)

Also writes plots/compare_base/base_chisq_summary.csv with IDL and XRTpy
chi-square side by side, plus which solver won and whether IDL stalled.

Pure load-and-plot -- no solving.

Usage (from synthetic_dem_validation/):
    python compare_synthetic_base_3way.py --dem 0     # one DEM
    python compare_synthetic_base_3way.py             # all found
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
PLOT_DIR = BASE_DIR / "plots" / "compare_base"

COLOR_TRUE = "red"
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"
STALL_THRESHOLD = 100.0  # chi2 above this = flagged as an IDL stall


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


# ── Loaders ───────────────────────────────────────────────────────────────────
def load_idl_base(dem_id):
    """Returns (logT, dem_base, chisq, filters, obs_val, obs_err)."""
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
    """
    Returns (logT, dem_base, chisq). Prefers the MC1000 npz (row 0 = base);
    falls back to a base-only npz if that's what exists.
    """
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
def plot_one_dem(dem_id):
    idl_logT, idl_base, idl_chisq, filters, obs_val, obs_err = load_idl_base(dem_id)
    py_logT, py_base, py_chisq = load_xrtpy_base(dem_id)

    if not np.allclose(idl_logT, py_logT, atol=1e-6):
        raise ValueError(f"IDL and XRTpy logT grids differ for DEM {dem_id}!")
    logT = py_logT

    fig, ax = plt.subplots(figsize=(11, 7))

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
    ax.set_ylim(10, 26)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.87])
    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    out_png = PLOT_DIR / f"compare_base_dem_idx{dem_id}_allfilters.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return idl_chisq, py_chisq


# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(p: Path) -> int:
    m = re.search(r"idx(\d+)_", p.name)
    return int(m.group(1)) if m else -1


def main():
    parser = argparse.ArgumentParser(description="3-way base DEM comparison")
    parser.add_argument("--dem", type=int, default=None,
                        help="Single DEM index (e.g., 0). Omit for all found.")
    args = parser.parse_args()

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

    print(f"Comparing {len(dem_ids)} DEM(s).")
    results = []
    for k, dem_id in enumerate(dem_ids, 1):
        idl_chi, py_chi = plot_one_dem(dem_id)
        results.append((dem_id, idl_chi, py_chi))
        winner = "XRTpy" if py_chi < idl_chi else "IDL"
        print(f"[{k:3d}/{len(dem_ids)}] DEM {dem_id:3d}  "
              f"IDL {idl_chi:12.4f}   XRTpy {py_chi:10.4f}   -> {winner}")

    if len(results) > 1:
        csv_path = PLOT_DIR / "base_chisq_summary.csv"
        with open(csv_path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["dem_id", "idl_base_chisq", "xrtpy_base_chisq",
                        "winner", "idl_stalled"])
            for dem_id, idl_chi, py_chi in sorted(results):
                winner = "XRTpy" if py_chi < idl_chi else "IDL"
                stalled = "yes" if idl_chi > STALL_THRESHOLD else "no"
                w.writerow([dem_id, f"{idl_chi:.6f}", f"{py_chi:.6f}",
                            winner, stalled])
        print(f"\nSummary saved: {csv_path}")

        idl = np.array([r[1] for r in results])
        py = np.array([r[2] for r in results])
        n_xrtpy_wins = int(np.sum(py < idl))
        n_idl_stall = int(np.sum(idl > STALL_THRESHOLD))
        print(f"XRTpy lower chi2 in {n_xrtpy_wins}/{len(results)} DEMs")
        print(f"IDL stalls (chi2 > {STALL_THRESHOLD:.0f}): "
              f"{n_idl_stall}/{len(results)} DEMs")
        print(f"IDL   chi2 median {np.median(idl):.2f}   max {idl.max():.2f}")
        print(f"XRTpy chi2 median {np.median(py):.2f}   max {py.max():.2f}")

    print("\nDone.")


if __name__ == "__main__":
    main()

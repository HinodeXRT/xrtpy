"""
Plot IDL DEM + 1000 Monte Carlo vs synthetic truth for the 162-DEM campaign
(IDL only; no XRTpy).

For each DEM index, loads:
    - True DEM:  data/synthetic_dems_data/DEM_{id}.txt
    - IDL base:  data/all_IDL_synthetic_dem_data/
                     idl_synthetic_dem_idx{id}_allfilters_base.sav
    - IDL MC:    data/all_IDL_MC_synthetic_dem_data/
                     idl_synthetic_dem_idx{id}_allfilters_MC1000.sav

and produces one plot per DEM:
    - IDL MC ensemble : all 1000 curves, faint orange
    - True DEM        : red solid
    - IDL base        : orange solid (chi-square in legend)

The base curve comes from the separate base .sav, NOT column 0 of the MC
output. Filter names / intensities / errors are read from the base .sav
(obs_index, obs_val, obs_err) and listed under the title.

Saved to plots/idl_mc1000/idl_synthetic_dem_idx{id}_allfilters_MC1000.png

Pure load-and-plot -- no solving.

Usage (from synthetic_dem_validation/):
    python plot_idl_synthetic_mc1000.py --dem 0     # one DEM
    python plot_idl_synthetic_mc1000.py             # all found
"""

import argparse
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
IDL_MC_DIR = BASE_DIR / "data" / "all_IDL_MC_synthetic_dem_data"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
PLOT_DIR = BASE_DIR / "plots" / "idl_mc1000"

N_MC = 1000
COLOR_TRUE = "red"
COLOR_IDL = "orange"


# ── Helpers ───────────────────────────────────────────────────────────────────
def _log10_dem(dem, floor=1e-40):
    return np.log10(np.maximum(np.asarray(dem, dtype=float), floor))


def _runs_by_T(dem, nT):
    """Return DEM array shaped (n_runs, nT) regardless of IDL orientation."""
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


def load_idl_mc(dem_id, logT):
    """MC DEM shaped (n_runs, nT), dropping column 0 (MC base column)."""
    path = IDL_MC_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_MC{N_MC}.sav"
    data = readsav(str(path), python_dict=True)
    dem_all = _runs_by_T(data["dem_out"], len(logT))
    return dem_all[1:]


def load_true_dem(dem_id):
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    if not path.exists():
        return None
    d = np.loadtxt(path)
    return d[:, 0], _log10_dem(d[:, 1])


# ── Plot ──────────────────────────────────────────────────────────────────────
def plot_one_dem(dem_id):
    logT, idl_base, idl_chisq, filters, obs_val, obs_err = load_idl_base(dem_id)
    idl_mc = load_idl_mc(dem_id, logT)

    fig, ax = plt.subplots(figsize=(11, 7))

    # All 1000 MC curves, faint orange
    for row in idl_mc:
        ax.step(logT, _log10_dem(row), where="mid",
                color=COLOR_IDL, alpha=0.02, linewidth=0.8)

    # True DEM
    true_dem = load_true_dem(dem_id)
    if true_dem is not None:
        ax.step(true_dem[0], true_dem[1], where="mid",
                color=COLOR_TRUE, linewidth=2.5, label="Synthetic DEM")

    # IDL base on top
    ax.step(logT, _log10_dem(idl_base), where="mid",
            color=COLOR_IDL, linewidth=2.5,
            label=f"IDL base  |  χ² = {idl_chisq:.4f}")

    ax.set_title(
        f"IDL DEM + {idl_mc.shape[0]} MC vs Synthetic — DEM Index {dem_id} "
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
    out_png = PLOT_DIR / f"idl_synthetic_dem_idx{dem_id}_allfilters_MC{N_MC}.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)
    return idl_chisq


# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(p: Path) -> int:
    m = re.search(r"idx(\d+)_", p.name)
    return int(m.group(1)) if m else -1


def main():
    parser = argparse.ArgumentParser(description="IDL MC1000 vs synthetic plots")
    parser.add_argument("--dem", type=int, default=None,
                        help="Single DEM index (e.g., 0). Omit for all found.")
    args = parser.parse_args()

    mc_files = sorted(
        IDL_MC_DIR.glob(f"idl_synthetic_dem_idx*_allfilters_MC{N_MC}.sav"),
        key=_dem_id_from_path,
    )
    if not mc_files:
        raise FileNotFoundError(f"No IDL MC .sav files found in {IDL_MC_DIR}")

    dem_ids = [_dem_id_from_path(p) for p in mc_files]
    if args.dem is not None:
        if args.dem not in dem_ids:
            raise FileNotFoundError(f"No IDL MC .sav for DEM {args.dem}")
        dem_ids = [args.dem]

    print(f"Plotting {len(dem_ids)} DEM(s).")
    for k, dem_id in enumerate(dem_ids, 1):
        chisq = plot_one_dem(dem_id)
        print(f"[{k:3d}/{len(dem_ids)}] DEM {dem_id:3d}  base chi2 = {chisq:12.4f}")

    print("\nDone.")


if __name__ == "__main__":
    main()

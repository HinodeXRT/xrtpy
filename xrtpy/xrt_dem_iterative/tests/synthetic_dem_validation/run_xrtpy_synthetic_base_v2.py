"""
XRTpy BASE DEM runner (v2) for the corrected-intensity 162-DEM campaign.

Reads intensity files from:
    data/synthetic_initial_condition_dems_data_with_bands/
(single source folder; C-poly/Al-thick skipped at read time)

For each DEM:
    - Run XRTDEMIterative, monte_carlo_runs=0 (base solution only)
    - Save xrtpy_output/xrtpy_dem_idx{id}_allfilters_base.npz
    - Save plot: true DEM (red) + XRTpy base (blue dashed, chi-square in legend)

At the end (full run only), writes a tidy chi-square CSV:
    chisq_csv/xrtpy_base_chisq.csv   with columns: dem_id, chisq
rebuilt from every base npz on disk, so it is complete even across resumes.

RESUME: a DEM whose base .npz exists is skipped (--overwrite to force).

Usage (from synthetic_dem_validation/):
    python run_xrtpy_synthetic_base_v2.py --dem 0     # one DEM (start here)
    python run_xrtpy_synthetic_base_v2.py             # all 162
    python run_xrtpy_synthetic_base_v2.py --overwrite # redo everything
"""

import argparse
import csv
import re
import time
import warnings
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from xrtpy.response.tools import generate_temperature_responses
from xrtpy.xrt_dem_iterative import XRTDEMIterative

# ── Configuration ─────────────────────────────────────────────────────────────
BASE_DIR = Path(__file__).parent
INTENSITY_DIR = BASE_DIR / "data" / "synthetic_initial_condition_dems_data_with_bands"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
OUT_DIR = BASE_DIR / "xrtpy_output"
PLOT_DIR = BASE_DIR / "plots"
CHISQ_DIR = BASE_DIR / "chisq_csv"

OBSERVATION_DATE = "2008-10-10T00:00:00"
EXCLUDED_FILTERS = {"C-poly/Al-thick"}

COLOR_TRUE = "red"
COLOR_XRTPY = "#1E90FF"


# ── Reader ────────────────────────────────────────────────────────────────────
def read_intensity_file(txt_path: Path):
    """Read intensity file; skip excluded filters AND the >2% band line."""
    lines = txt_path.read_text().strip().splitlines()
    dem_id = int(lines[0].split(":")[1].strip())

    filters, intensities, uncertainties = [], [], []
    for line in lines:
        stripped = line.strip()
        low = stripped.lower()
        # Skip headers, units, the contribution line, the column header,
        # and any continuation lines of the contribution list.
        if (low.startswith("dem index") or low.startswith("units")
                or low.startswith("temperature") or low.startswith("filter")
                or "-" in stripped[:0]):  # placeholder, no-op
            continue
        parts = stripped.split()
        if len(parts) != 4:
            continue  # skip date line, band continuation lines, blanks
        name = parts[0]
        # A valid filter row has a non-numeric first token
        try:
            float(name)
            continue  # first token is a number -> band continuation line
        except ValueError:
            pass
        if name in EXCLUDED_FILTERS:
            continue
        try:
            i_val = float(parts[1]); e_val = float(parts[2]); float(parts[3])
        except ValueError:
            continue
        filters.append(name)
        intensities.append(i_val)
        uncertainties.append(e_val)

    return {
        "dem_id": dem_id,
        "filters": filters,
        "intensities": np.array(intensities, dtype=float),
        "uncertainties": np.array(uncertainties, dtype=float),
    }


def read_true_dem(dem_id: int):
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"
    if not path.exists():
        return None
    d = np.loadtxt(path)
    return d[:, 0], np.log10(np.clip(d[:, 1], 1e-40, None))


# ── Runner ────────────────────────────────────────────────────────────────────
def run_one_dem(txt_path: Path, responses_cache: dict):
    case = read_intensity_file(txt_path)
    dem_id = case["dem_id"]
    filters = case["filters"]

    key = tuple(filters)
    if key not in responses_cache:
        responses_cache[key] = generate_temperature_responses(
            filters, OBSERVATION_DATE
        )
    responses = responses_cache[key]

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=UserWarning)
        solver = XRTDEMIterative(
            observed_channel=filters,
            observed_intensities=case["intensities"],
            temperature_responses=responses,
            intensity_uncertainties=case["uncertainties"],
            monte_carlo_runs=0,
        )
        solver.solve()

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    out_npz = OUT_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
    np.savez(
        out_npz,
        dem_id=dem_id,
        filters=np.array(filters),
        observation_date=OBSERVATION_DATE,
        intensities=case["intensities"],
        uncertainties=case["uncertainties"],
        logT=solver.logT,
        dem=solver.dem,
        chisq=solver.chisq,
        modeled_intensities=solver.modeled_intensities,
    )

    PLOT_DIR.mkdir(parents=True, exist_ok=True)
    plot_one_dem(solver, dem_id, filters, case["intensities"],
                 case["uncertainties"])
    return dem_id, float(solver.chisq)


def plot_one_dem(solver, dem_id, filters, intensities, uncertainties):
    logT = solver.logT
    log_dem = np.log10(np.clip(solver.dem, 1e-40, None))

    fig, ax = plt.subplots(figsize=(11, 7))
    true_dem = read_true_dem(dem_id)
    if true_dem is not None:
        ax.step(true_dem[0], true_dem[1], where="mid",
                color=COLOR_TRUE, linewidth=2.5, label="Synthetic DEM")
    ax.step(logT, log_dem, where="mid",
            color=COLOR_XRTPY, linewidth=2.5, linestyle="--",
            label=f"XRTpy base  |  χ² = {solver.chisq:.4f}")

    ax.set_title(
        f"XRTpy Base DEM vs Synthetic — DEM Index {dem_id} "
        f"(all {len(filters)} filters)",
        fontsize=13, pad=34,
    )
    entries = [
        f"{f} [{i:.4g} ± {e:.3g}]"
        for f, i, e in zip(filters, intensities, uncertainties, strict=True)
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
    ax.set_ylim(14, 26)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.87])
    out_png = PLOT_DIR / f"xrtpy_synthetic_dem_idx{dem_id}_allfilters_base.png"
    fig.savefig(out_png, dpi=150)
    plt.close(fig)


def write_base_chisq_csv():
    """Rebuild xrtpy_base_chisq.csv from all base npz on disk."""
    npz_files = sorted(
        OUT_DIR.glob("xrtpy_dem_idx*_allfilters_base.npz"),
        key=lambda p: int(re.search(r"idx(\d+)_", p.name).group(1)),
    )
    if not npz_files:
        return
    CHISQ_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = CHISQ_DIR / "xrtpy_base_chisq.csv"
    with open(csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["dem_id", "chisq"])
        for p in npz_files:
            d = np.load(p, allow_pickle=True)
            w.writerow([int(d["dem_id"]), f"{float(d['chisq']):.6f}"])
    print(f"\nChi-square CSV ({len(npz_files)} DEMs): {csv_path}")


# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(p: Path) -> int:
    m = re.search(r"DEM_(\d+)\.txt$", p.name)
    return int(m.group(1)) if m else -1


def main():
    parser = argparse.ArgumentParser(description="XRTpy base runner v2")
    parser.add_argument("--dem", type=int, default=None,
                        help="Single DEM index (e.g., 0). Omit for all.")
    parser.add_argument("--overwrite", action="store_true",
                        help="Re-solve even if the base .npz exists.")
    args = parser.parse_args()

    txt_files = sorted(
        INTENSITY_DIR.glob("XRT_intensities_DEM_*.txt"),
        key=_dem_id_from_path,
    )
    if not txt_files:
        raise FileNotFoundError(f"No intensity files in {INTENSITY_DIR}")

    if args.dem is not None:
        txt_files = [p for p in txt_files if _dem_id_from_path(p) == args.dem]
        if not txt_files:
            raise FileNotFoundError(f"No file for DEM index {args.dem}")

    todo, skipped = [], 0
    for p in txt_files:
        dem_id = _dem_id_from_path(p)
        out_npz = OUT_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
        if out_npz.exists() and not args.overwrite:
            skipped += 1
        else:
            todo.append(p)

    print(f"Found {len(txt_files)} file(s); {skipped} done, {len(todo)} to solve.")
    if not todo:
        print("Nothing to do. Use --overwrite to re-solve.")
        # still (re)write the CSV so it reflects disk
        if args.dem is None:
            write_base_chisq_csv()
        return

    responses_cache = {}
    t0 = time.time()
    for k, p in enumerate(todo, 1):
        dem_id, chisq = run_one_dem(p, responses_cache)
        elapsed = time.time() - t0
        print(f"[{k:3d}/{len(todo)}] DEM {dem_id:3d}  chi2 = {chisq:12.4f}   "
              f"({elapsed/60:.1f} min)", flush=True)

    # Write the tidy chi-square CSV for the full population
    if args.dem is None:
        write_base_chisq_csv()

    print("\nDone.")


if __name__ == "__main__":
    main()
"""
Compare synthetic truth, IDL base DEM results, and XRTpy base DEM results
with contribution-band shading.

For each synthetic DEM case, the plot contains:

- synthetic truth DEM;
- IDL base DEM;
- XRTpy base DEM;
- shaded temperature ranges corresponding to the contribution-band data.

Contiguous qualifying temperature bins are grouped into separate shaded
regions so multi-component DEMs retain gaps between constrained regions.

Inputs
------
Contribution bands:
    data/contribution_bands.csv

Synthetic truth:
    data/synthetic_dems_data/DEM_{id}.txt

IDL base results:
    data/all_IDL_synthetic_dem_data/
    idl_synthetic_dem_idx{id}_allfilters_base.sav

XRTpy base results:
    data/output_data/xrtpy_base_dem_results/
    xrtpy_dem_idx{id}_allfilters_base.npz

Outputs
-------
    plots/idl_xrtpy_synthetic_base_dem_comparison_with_contribution_bands/
    compare_base_dem_idx{id}_allfilters.png

Usage from ``synthetic_dem_validation/``:

    python compare_base_3way_banded.py --dem 137
    python compare_base_3way_banded.py
"""

import argparse
import csv
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.io import readsav


BASE_DIR = Path(__file__).parent

IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"

XRTPY_BASE_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "xrtpy_base_dem_results"
)

BANDS_CSV = BASE_DIR / "data" / "contribution_bands.csv"

PLOT_DIR = (
    BASE_DIR
    / "plots"
    / "idl_xrtpy_synthetic_base_dem_comparison_with_contribution_bands"
)


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
    """Load the XRTpy base DEM result for one synthetic case."""
    path = (
        XRTPY_BASE_DIR
        / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"No XRTpy base NPZ for DEM {dem_id}: {path}"
        )

    with np.load(path) as result:
        log_temperature = np.asarray(result["logT"], dtype=float)
        dem = np.asarray(result["dem"], dtype=float)
        chisq = float(result["chisq"])

    return log_temperature, dem, chisq

def load_true_dem(dem_id):
    """Load the synthetic truth DEM for one case."""
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"

    if not path.exists():
        return None

    data = np.loadtxt(path)

    log_temperature = data[:, 0]
    log_dem = _log10_dem(data[:, 1])

    return log_temperature, log_dem



# ── Plot ──────────────────────────────────────────────────────────────────────
def plot_one_dem(dem_id, bands):
    """Plot truth, IDL base, and XRTpy base DEMs with contribution bands."""
    (
        idl_log_temperature,
        idl_base_dem,
        idl_chisq,
        filters,
        observed_values,
        observed_errors,
    ) = load_idl_base(dem_id)

    (
        xrtpy_log_temperature,
        xrtpy_base_dem,
        xrtpy_chisq,
    ) = load_xrtpy_base(dem_id)

    if not np.allclose(
        idl_log_temperature,
        xrtpy_log_temperature,
        atol=1e-6,
    ):
        raise ValueError(
            f"IDL and XRTpy logT grids differ for DEM {dem_id}."
        )

    log_temperature = xrtpy_log_temperature

    fig, ax = plt.subplots(figsize=(11, 7))

    spans = contiguous_spans(
        bands.get(dem_id, [])
    )

    for span_index, (lower, upper) in enumerate(spans):
        ax.axvspan(
            lower,
            upper,
            color=BAND_COLOR,
            alpha=BAND_ALPHA,
            label=">2% contribution" if span_index == 0 else None,
        )

    true_dem = load_true_dem(dem_id)

    if true_dem is not None:
        true_log_temperature, true_log_dem = true_dem

        ax.step(
            true_log_temperature,
            true_log_dem,
            where="mid",
            color=COLOR_TRUE,
            linewidth=2.5,
            label="Synthetic DEM",
        )

    ax.step(
        log_temperature,
        _log10_dem(idl_base_dem),
        where="mid",
        color=COLOR_IDL,
        linewidth=2.5,
        label=f"IDL base  |  χ² = {idl_chisq:.4f}",
    )

    ax.step(
        log_temperature,
        _log10_dem(xrtpy_base_dem),
        where="mid",
        color=COLOR_XRTPY,
        linewidth=2.5,
        linestyle="--",
        label=f"XRTpy base  |  χ² = {xrtpy_chisq:.4f}",
    )

    ax.set_title(
        f"Base DEM Comparison — IDL vs XRTpy vs Synthetic — "
        f"DEM Index {dem_id} "
        f"(all {len(filters)} filters)",
        fontsize=13,
        pad=34,
    )

    filter_entries = [
        f"{filter_name} [{value:.4g} ± {error:.3g}]"
        for filter_name, value, error in zip(
            filters,
            observed_values,
            observed_errors,
            strict=True,
        )
    ]

    filters_per_line = 4
    filter_label = "\n".join(
        ",   ".join(
            filter_entries[start : start + filters_per_line]
        )
        for start in range(0, len(filter_entries), filters_per_line)
    )

    fig.text(
        0.5,
        0.945,
        filter_label,
        ha="center",
        va="top",
        fontsize=7.5,
    )

    ax.set_xlabel(r"log$_{10}$ T  [K]", fontsize=12)
    ax.set_ylabel(r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]", fontsize=12)
    ax.set_xlim(log_temperature.min(), log_temperature.max())
    ax.set_ylim(18, 26)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.87])

    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    out_png = (
        PLOT_DIR
        / f"compare_base_dem_idx{dem_id}_allfilters.png"
    )

    fig.savefig(out_png, dpi=150)
    plt.close(fig)

    return idl_chisq, xrtpy_chisq, len(spans)



# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(path: Path) -> int:
    """Extract the DEM index from an IDL result filename."""
    match = re.search(r"idx(\d+)_", path.name)

    if match is None:
        raise ValueError(
            f"Could not extract DEM index from filename: {path.name}"
        )

    return int(match.group(1))


def main():
    """Generate contribution-banded base DEM comparison plots."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare synthetic truth, IDL base, and XRTpy base DEMs "
            "with contribution-band shading."
        )
    )
    parser.add_argument(
        "--dem",
        type=int,
        default=None,
        help="Plot a single DEM index. Omit to process all available cases.",
    )
    args = parser.parse_args()

    bands = load_bands()

    sav_files = sorted(
        IDL_BASE_DIR.glob(
            "idl_synthetic_dem_idx*_allfilters_base.sav"
        ),
        key=_dem_id_from_path,
    )

    if not sav_files:
        raise FileNotFoundError(
            f"No IDL base SAV files found in {IDL_BASE_DIR}"
        )

    dem_ids = [
        _dem_id_from_path(path)
        for path in sav_files
    ]

    if args.dem is not None:
        if args.dem not in dem_ids:
            raise FileNotFoundError(
                f"No IDL base SAV file found for DEM {args.dem}"
            )

        dem_ids = [args.dem]

    n_no_band = 0

    print(f"Plotting {len(dem_ids)} DEM(s).")

    for index, dem_id in enumerate(dem_ids, start=1):
        idl_chisq, xrtpy_chisq, n_spans = plot_one_dem(
            dem_id,
            bands,
        )

        if n_spans == 0:
            n_no_band += 1

        print(
            f"[{index:3d}/{len(dem_ids)}] "
            f"DEM {dem_id:3d}  "
            f"IDL {idl_chisq:12.4f}   "
            f"XRTpy {xrtpy_chisq:10.4f}   "
            f"bands: {n_spans}"
        )

    if n_no_band:
        print(
            f"\nNote: {n_no_band} DEM(s) plotted without shading "
            f"(not in {BANDS_CSV.name})."
        )

    print("\nDone.")


if __name__ == "__main__":
    main()

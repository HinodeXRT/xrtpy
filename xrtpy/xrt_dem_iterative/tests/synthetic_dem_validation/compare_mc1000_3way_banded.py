"""
Compare synthetic truth, IDL results, and XRTpy results with Monte Carlo
ensembles and contribution-band shading.

For each synthetic DEM case, the figure contains:

- contribution-band regions;
- IDL Monte Carlo ensemble;
- XRTpy Monte Carlo ensemble;
- synthetic truth DEM;
- standalone IDL base DEM;
- standalone XRTpy base DEM.

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

IDL Monte Carlo results:
    data/all_IDL_MC_synthetic_dem_data/
    idl_synthetic_dem_idx{id}_allfilters_MC1000.sav

XRTpy base results:
    data/output_data/xrtpy_base_dem_results/
    xrtpy_dem_idx{id}_allfilters_base.npz

XRTpy Monte Carlo results:
    data/output_data/xrtpy_monte_carlo_dem_results/
    xrtpy_dem_idx{id}_allfilters_MC1000.npz

Outputs
-------
    plots/idl_xrtpy_synthetic_monte_carlo_1000_run_comparison_with_contribution_bands/
    compare_mc1000_dem_idx{id}_allfilters.png

Usage from ``synthetic_dem_validation/``:

    python compare_mc1000_3way_banded.py --dem 137
    python compare_mc1000_3way_banded.py
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
IDL_MC_DIR = BASE_DIR / "data" / "all_IDL_MC_synthetic_dem_data"
TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
XRTPY_DIR = BASE_DIR / "data" / "output_data" / "xrtpy_monte_carlo_dem_results"
BANDS_CSV = BASE_DIR / "data" / "contribution_bands.csv"
PLOT_DIR = BASE_DIR / "plots" / "idl_xrtpy_synthetic_monte_carlo_1000_run_comparison_with_contribution_bands"

XRTPY_BASE_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "xrtpy_base_dem_results"
)

XRTPY_MC_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "xrtpy_monte_carlo_dem_results"
)

N_MC = 1000
COLOR_TRUE = "red"
COLOR_IDL = "orange"
COLOR_XRTPY = "#1E90FF"
ALPHA_MC = 0.02
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
    """Load contribution-band temperatures keyed by DEM index."""
    bands = {}

    if not BANDS_CSV.exists():
        print(
            f"Note: {BANDS_CSV} not found -- "
            "plotting without contribution-band shading."
        )
        return bands

    with BANDS_CSV.open(newline="") as csv_file:
        reader = csv.DictReader(csv_file)

        for row in reader:
            temps = (
                [float(value) for value in row["temps"].split()]
                if row["temps"]
                else []
            )

            bands[int(row["dem_id"])] = temps

    return bands


def contiguous_spans(
    temps,
    step=GRID_STEP,
):
    """Group contribution temperatures into contiguous logT spans."""
    if not temps:
        return []

    temperatures = np.sort(
        np.asarray(temps, dtype=float)
    )

    spans = []

    start = temperatures[0]
    previous = temperatures[0]

    for temperature in temperatures[1:]:
        if temperature - previous > step * 1.5:
            spans.append(
                (start, previous)
            )
            start = temperature

        previous = temperature

    spans.append(
        (start, previous)
    )

    half_step = step / 2.0

    return [
        (
            lower - half_step,
            upper + half_step,
        )
        for lower, upper in spans
    ]


# ── Loaders ───────────────────────────────────────────────────────────────────
def load_idl_base(dem_id):
    """Load the standalone IDL base DEM result."""
    path = (
        IDL_BASE_DIR
        / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"
    )

    data = readsav(
        str(path),
        python_dict=True,
    )

    log_temperature = np.asarray(
        data["logt_out"],
        dtype=float,
    ).ravel()

    base_dem = _runs_by_T(
        data["dem_out"],
        len(log_temperature),
    )[0]

    chisq = float(
        np.atleast_1d(data["chisq"])[0]
    )

    filters = _decode(
        data["obs_index"]
    )

    observed_values = np.atleast_1d(
        data["obs_val"]
    ).astype(float)

    observed_errors = np.atleast_1d(
        data["obs_err"]
    ).astype(float)

    return (
        log_temperature,
        base_dem,
        chisq,
        filters,
        observed_values,
        observed_errors,
    )


def load_idl_mc(
    dem_id,
    log_temperature,
):
    """Load the IDL Monte Carlo DEM realizations, excluding row 0."""
    path = (
        IDL_MC_DIR
        / f"idl_synthetic_dem_idx{dem_id}_allfilters_MC{N_MC}.sav"
    )

    data = readsav(
        str(path),
        python_dict=True,
    )

    all_dem_runs = _runs_by_T(
        data["dem_out"],
        len(log_temperature),
    )

    return all_dem_runs[1:]


def load_xrtpy_base(dem_id):
    """Load the standalone XRTpy base DEM result."""
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
        base_dem = np.asarray(result["dem"], dtype=float)
        chisq = float(result["chisq"])

    return log_temperature, base_dem, chisq


def load_xrtpy_mc(dem_id):
    """Load the XRTpy Monte Carlo DEM realizations."""
    path = (
        XRTPY_MC_DIR
        / f"xrtpy_dem_idx{dem_id}_allfilters_MC{N_MC}.npz"
    )

    if not path.exists():
        raise FileNotFoundError(
            f"No XRTpy MC NPZ for DEM {dem_id}: {path}"
        )

    with np.load(path) as result:
        log_temperature = np.asarray(result["logT"], dtype=float)
        mc_dem = np.asarray(result["mc_dem"], dtype=float)

    return log_temperature, mc_dem[1:]



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
    """Plot truth, base DEMs, MC ensembles, and contribution bands."""
    (
        idl_log_temperature,
        idl_base_dem,
        idl_chisq,
        filters,
        observed_values,
        observed_errors,
    ) = load_idl_base(dem_id)

    idl_mc = load_idl_mc(
        dem_id,
        idl_log_temperature,
    )

    (
        xrtpy_base_log_temperature,
        xrtpy_base_dem,
        xrtpy_chisq,
    ) = load_xrtpy_base(dem_id)

    (
        xrtpy_mc_log_temperature,
        xrtpy_mc,
    ) = load_xrtpy_mc(dem_id)

    if not np.allclose(
        idl_log_temperature,
        xrtpy_base_log_temperature,
        atol=1e-6,
    ):
        raise ValueError(
            f"IDL and XRTpy base logT grids differ for DEM {dem_id}."
        )

    if not np.allclose(
        idl_log_temperature,
        xrtpy_mc_log_temperature,
        atol=1e-6,
    ):
        raise ValueError(
            f"IDL and XRTpy MC logT grids differ for DEM {dem_id}."
        )

    log_temperature = xrtpy_base_log_temperature

    fig, ax = plt.subplots(figsize=(11, 7))

    # Contribution bands are drawn first so they remain behind all DEM curves.
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

    # Monte Carlo ensembles.
    for mc_dem in idl_mc:
        ax.step(
            log_temperature,
            _log10_dem(mc_dem),
            where="mid",
            color=COLOR_IDL,
            alpha=ALPHA_MC,
            linewidth=0.8,
        )

    for mc_dem in xrtpy_mc:
        ax.step(
            log_temperature,
            _log10_dem(mc_dem),
            where="mid",
            color=COLOR_XRTPY,
            alpha=ALPHA_MC,
            linewidth=0.8,
        )

    # Synthetic truth.
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

    # Standalone base solutions are plotted on top of the MC ensembles.
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
        f"DEM + {N_MC} MC — IDL vs XRTpy vs Synthetic — "
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
        for start in range(
            0,
            len(filter_entries),
            filters_per_line,
        )
    )

    fig.text(
        0.5,
        0.945,
        filter_label,
        ha="center",
        va="top",
        fontsize=7.5,
    )

    ax.set_xlabel(
        r"log$_{10}$ T  [K]",
        fontsize=12,
    )
    ax.set_ylabel(
        r"log$_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]",
        fontsize=12,
    )

    ax.set_xlim(
        log_temperature.min(),
        log_temperature.max(),
    )
    ax.set_ylim(16, 26)

    ax.grid(
        visible=True,
        alpha=0.3,
    )
    ax.legend(
        fontsize=11,
        loc="best",
    )

    fig.tight_layout(
        rect=[0, 0, 1, 0.87]
    )

    PLOT_DIR.mkdir(
        parents=True,
        exist_ok=True,
    )

    out_png = (
        PLOT_DIR
        / f"compare_mc{N_MC}_dem_idx{dem_id}_allfilters.png"
    )

    fig.savefig(
        out_png,
        dpi=150,
    )
    plt.close(fig)

    return (
        idl_chisq,
        xrtpy_chisq,
        len(spans),
    )



# ── Main ──────────────────────────────────────────────────────────────────────

def _dem_id_from_path(path: Path) -> int:
    """Extract the DEM index from a result filename."""
    match = re.search(r"idx(\d+)_", path.name)

    if match is None:
        raise ValueError(
            f"Could not extract DEM index from filename: {path.name}"
        )

    return int(match.group(1))


def main():
    """Generate MC1000 comparison plots with contribution-band shading."""
    parser = argparse.ArgumentParser(
        description=(
            "Compare synthetic truth, IDL, and XRTpy base and "
            "Monte Carlo DEM results with contribution-band shading."
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

    idl_ids = {
        _dem_id_from_path(path)
        for path in IDL_MC_DIR.glob(
            f"idl_synthetic_dem_idx*_allfilters_MC{N_MC}.sav"
        )
    }

    xrtpy_ids = {
        _dem_id_from_path(path)
        for path in XRTPY_MC_DIR.glob(
            f"xrtpy_dem_idx*_allfilters_MC{N_MC}.npz"
        )
    }

    dem_ids = sorted(
        idl_ids & xrtpy_ids
    )

    if not dem_ids:
        raise FileNotFoundError(
            "No DEMs with both IDL and XRTpy MC results were found."
        )

    missing = sorted(
        (idl_ids | xrtpy_ids)
        - (idl_ids & xrtpy_ids)
    )

    if missing:
        preview = missing[:15]
        suffix = " ..." if len(missing) > 15 else ""

        print(
            f"Note: {len(missing)} DEM(s) skipped "
            f"(one solver only): {preview}{suffix}"
        )

    if args.dem is not None:
        if args.dem not in dem_ids:
            raise FileNotFoundError(
                f"DEM {args.dem} is missing IDL and/or XRTpy MC output."
            )

        dem_ids = [args.dem]

    print(f"Comparing {len(dem_ids)} DEM(s).")

    for index, dem_id in enumerate(dem_ids, start=1):
        idl_chisq, xrtpy_chisq, n_spans = plot_one_dem(
            dem_id,
            bands,
        )

        print(
            f"[{index:3d}/{len(dem_ids)}] "
            f"DEM {dem_id:3d}  "
            f"IDL {idl_chisq:12.4f}   "
            f"XRTpy {xrtpy_chisq:10.4f}   "
            f"bands: {n_spans}"
        )

    print("\nDone.")

if __name__ == "__main__":
    main()
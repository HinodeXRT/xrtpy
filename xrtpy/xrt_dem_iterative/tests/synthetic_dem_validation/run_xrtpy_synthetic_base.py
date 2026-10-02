"""
Run XRTpy base DEM inversions for the 162-case synthetic validation campaign.

Synthetic XRT intensities are read from:

    data/synthetic_initial_condition_dems_data_with_bands/

The C-poly/Al-thick channel is excluded when the intensity files are read.

For each synthetic DEM case, this script:

- generates the XRT temperature responses for the selected filters;
- runs ``XRTDEMIterative`` with ``monte_carlo_runs=0``;
- saves the base inversion result to:

      data/output_data/xrtpy_base_dem_results/

- saves a comparison plot of the synthetic truth and XRTpy base solution to:

      plots/xrtpy_base_dem_vs_synthetic_truth/

For a full campaign run, the script also rebuilds:

    data/output_data/chi_square_comparison_tables/xrtpy_base_chisq.csv

The chi-square table is generated from all saved base NPZ files on disk, so it
remains complete when a campaign is resumed.

Existing base NPZ files are skipped unless ``--overwrite`` is supplied.

Usage from ``synthetic_dem_validation/``:

    python run_xrtpy_synthetic_base.py --dem 0
    python run_xrtpy_synthetic_base.py
    python run_xrtpy_synthetic_base.py --overwrite
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
OUT_DIR = BASE_DIR / "data" / "output_data" / "xrtpy_base_dem_results"
PLOT_DIR = BASE_DIR / "plots" / "xrtpy_base_dem_vs_synthetic_truth"
CHISQ_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "chi_square_comparison_tables"
)

OBSERVATION_DATE = "2008-10-10T00:00:00"
EXCLUDED_FILTERS = {"C-poly/Al-thick"}

COLOR_TRUE = "red"
COLOR_XRTPY = "#1E90FF"


# ── Reader ────────────────────────────────────────────────────────────────────
def read_intensity_file(txt_path: Path):
    """Read a synthetic XRT intensity file.

    The file contains metadata, contribution-band information, and filter
    intensity rows. Only valid filter rows are returned. Filters listed in
    ``EXCLUDED_FILTERS`` are skipped.

    Parameters
    ----------
    txt_path : pathlib.Path
        Path to the synthetic intensity text file.

    Returns
    -------
    dict
        Dictionary containing the DEM index, filter names, intensities, and
        intensity uncertainties.
    """
    lines = txt_path.read_text().strip().splitlines()
    dem_id = int(lines[0].split(":")[1].strip())

    filters = []
    intensities = []
    uncertainties = []

    for line in lines:
        stripped = line.strip()
        low = stripped.lower()

        # Skip metadata and table headers.
        if (
            low.startswith("dem index")
            or low.startswith("units")
            or low.startswith("temperature")
            or low.startswith("filter")
        ):
            continue

        parts = stripped.split()

        # Valid filter rows contain four columns.
        # This also skips blank lines, date lines, and contribution-band lines.
        if len(parts) != 4:
            continue

        name = parts[0]

        # Contribution-band continuation rows begin with a numeric value.
        try:
            float(name)
        except ValueError:
            pass
        else:
            continue

        if name in EXCLUDED_FILTERS:
            continue

        try:
            intensity = float(parts[1])
            uncertainty = float(parts[2])

            # Validate the fourth column even though it is not currently used.
            float(parts[3])
        except ValueError:
            continue

        filters.append(name)
        intensities.append(intensity)
        uncertainties.append(uncertainty)

    return {
        "dem_id": dem_id,
        "filters": filters,
        "intensities": np.asarray(intensities, dtype=float),
        "uncertainties": np.asarray(uncertainties, dtype=float),
    }


def read_true_dem(dem_id: int):
    """Load the synthetic truth DEM for a given DEM index.

    Parameters
    ----------
    dem_id : int
        Synthetic DEM index.

    Returns
    -------
    tuple[numpy.ndarray, numpy.ndarray] or None
        Temperature grid and log10 DEM values. Returns ``None`` if the
        corresponding truth DEM file does not exist.
    """
    path = TRUE_DEM_DIR / f"DEM_{dem_id}.txt"

    if not path.exists():
        return None

    data = np.loadtxt(path)

    log_temperature = data[:, 0]
    log_dem = np.log10(np.clip(data[:, 1], 1e-40, None))

    return log_temperature, log_dem

# ── Runner ────────────────────────────────────────────────────────────────────
def run_one_dem(txt_path: Path, responses_cache: dict):
    """Run the XRTpy base DEM inversion for one synthetic DEM case.

    Parameters
    ----------
    txt_path : pathlib.Path
        Path to the synthetic intensity file for one DEM case.
    responses_cache : dict
        Cache of temperature responses keyed by the ordered filter tuple.
        Cases using the same filter set reuse the previously generated
        responses.

    Returns
    -------
    tuple[int, float]
        DEM index and chi-square value for the base inversion.
    """
    case = read_intensity_file(txt_path)

    dem_id = case["dem_id"]
    filters = case["filters"]

    response_key = tuple(filters)

    if response_key not in responses_cache:
        responses_cache[response_key] = generate_temperature_responses(
            filters,
            OBSERVATION_DATE,
        )

    responses = responses_cache[response_key]

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

    plot_one_dem(
        solver,
        dem_id,
        filters,
        case["intensities"],
        case["uncertainties"],
    )

    return dem_id, float(solver.chisq)


def plot_one_dem(solver, dem_id, filters, intensities, uncertainties):
    """Plot the XRTpy base DEM against the known synthetic truth.

    Parameters
    ----------
    solver : XRTDEMIterative
        Solved XRTpy DEM inversion instance.
    dem_id : int
        Synthetic DEM index.
    filters : list[str]
        XRT filters used in the inversion.
    intensities : numpy.ndarray
        Observed synthetic intensities for the selected filters.
    uncertainties : numpy.ndarray
        Intensity uncertainties for the selected filters.
    """
    log_temperature = solver.logT
    log_dem = np.log10(np.clip(solver.dem, 1e-40, None))

    fig, ax = plt.subplots(figsize=(11, 7))

    true_dem = read_true_dem(dem_id)
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
        log_dem,
        where="mid",
        color=COLOR_XRTPY,
        linewidth=2.5,
        linestyle="--",
        label=f"XRTpy base  |  χ² = {solver.chisq:.4f}",
    )

    ax.set_title(
        f"XRTpy Base DEM vs Synthetic — DEM Index {dem_id} "
        f"(all {len(filters)} filters)",
        fontsize=13,
        pad=34,
    )

    filter_entries = [
        f"{filter_name} [{intensity:.4g} ± {uncertainty:.3g}]"
        for filter_name, intensity, uncertainty in zip(
            filters,
            intensities,
            uncertainties,
            strict=True,
        )
    ]

    filters_per_line = 4
    filter_label = "\n".join(
        ",   ".join(filter_entries[start : start + filters_per_line])
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
    ax.set_ylim(14, 26)
    ax.grid(visible=True, alpha=0.3)
    ax.legend(fontsize=11, loc="best")

    fig.tight_layout(rect=[0, 0, 1, 0.87])

    out_png = (
        PLOT_DIR
        / f"xrtpy_synthetic_dem_idx{dem_id}_allfilters_base.png"
    )
    fig.savefig(out_png, dpi=150)
    plt.close(fig)

def write_base_chisq_csv():
    """Rebuild the base chi-square summary CSV from saved NPZ results."""
    npz_files = sorted(
        OUT_DIR.glob("xrtpy_dem_idx*_allfilters_base.npz"),
        key=lambda path: int(re.search(r"idx(\d+)_", path.name).group(1)),
    )

    if not npz_files:
        return

    CHISQ_DIR.mkdir(parents=True, exist_ok=True)
    csv_path = CHISQ_DIR / "xrtpy_base_chisq.csv"

    with csv_path.open("w", newline="") as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(["dem_id", "chisq"])

        for npz_path in npz_files:
            with np.load(npz_path) as result:
                dem_id = int(result["dem_id"])
                chisq = float(result["chisq"])

            writer.writerow([dem_id, f"{chisq:.6f}"])

    print(f"\nChi-square CSV ({len(npz_files)} DEMs): {csv_path}")

# ── Main ──────────────────────────────────────────────────────────────────────
def _dem_id_from_path(path: Path) -> int:
    """Extract the synthetic DEM index from an intensity filename."""
    match = re.search(r"DEM_(\d+)\.txt$", path.name)

    if match is None:
        raise ValueError(f"Could not extract DEM index from filename: {path.name}")

    return int(match.group(1))


def main():
    """Run one or all synthetic base DEM inversions."""
    parser = argparse.ArgumentParser(
        description="Run XRTpy base inversions for the synthetic DEM campaign."
    )
    parser.add_argument(
        "--dem",
        type=int,
        default=None,
        help="Run a single DEM index. Omit to process all available DEMs.",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Re-run DEMs even when their base NPZ result already exists.",
    )
    args = parser.parse_args()

    txt_files = sorted(
        INTENSITY_DIR.glob("XRT_intensities_DEM_*.txt"),
        key=_dem_id_from_path,
    )

    if not txt_files:
        raise FileNotFoundError(
            f"No synthetic intensity files found in {INTENSITY_DIR}"
        )

    if args.dem is not None:
        txt_files = [
            path
            for path in txt_files
            if _dem_id_from_path(path) == args.dem
        ]

        if not txt_files:
            raise FileNotFoundError(
                f"No intensity file found for DEM index {args.dem}"
            )

    todo = []
    skipped = 0

    for txt_path in txt_files:
        dem_id = _dem_id_from_path(txt_path)
        out_npz = OUT_DIR / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"

        if out_npz.exists() and not args.overwrite:
            skipped += 1
        else:
            todo.append(txt_path)

    print(
        f"Found {len(txt_files)} file(s); "
        f"{skipped} done, {len(todo)} to solve."
    )

    if not todo:
        print("Nothing to do. Use --overwrite to re-solve.")

        if args.dem is None:
            write_base_chisq_csv()

        return

    responses_cache = {}
    start_time = time.time()

    for index, txt_path in enumerate(todo, start=1):
        dem_id, chisq = run_one_dem(txt_path, responses_cache)

        elapsed_minutes = (time.time() - start_time) / 60

        print(
            f"[{index:3d}/{len(todo)}] "
            f"DEM {dem_id:3d}  "
            f"chi2 = {chisq:12.4f}   "
            f"({elapsed_minutes:.1f} min)",
            flush=True,
        )

    if args.dem is None:
        write_base_chisq_csv()

    print("\nDone.")


if __name__ == "__main__":
    main()
"""
Compare IDL and XRTpy Monte Carlo DEM results for the synthetic validation
campaign.

Both solvers use the same row convention:

    row 0        = base, unperturbed solution
    rows 1..1000 = Monte Carlo realizations

For each DEM index available in both result sets, this script:

1. compares the IDL and XRTpy Monte Carlo DEM ensembles;
2. builds a long chi-square table containing every paired MC realization;
3. builds a per-DEM chi-square summary table;
4. generates per-DEM ensemble overlays;
5. generates campaign-level chi-square distribution plots.

Inputs
------
IDL Monte Carlo results:
    data/all_IDL_MC_synthetic_dem_data/
    idl_synthetic_dem_idx{N}_allfilters_MC1000.sav

XRTpy Monte Carlo results:
    data/output_data/xrtpy_monte_carlo_dem_results/
    xrtpy_dem_idx{N}_allfilters_MC1000.npz

Synthetic truth DEMs:
    data/synthetic_dems_data/
    DEM_{N}.txt

Outputs
-------
Chi-square tables:
    data/output_data/chi_square_comparison_tables/
    mc1000_chi_square_all_runs.csv
    mc1000_chi_square_summary.csv

Per-DEM overlays:
    plots/idl_xrtpy_synthetic_monte_carlo_1000_run_comparison/
    per_dem_monte_carlo_overlays/

Chi-square summary plots:
    plots/idl_xrtpy_synthetic_monte_carlo_1000_run_comparison/
    chi_square_summary_plots/
"""


from __future__ import annotations

from pathlib import Path
import re

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.io import readsav

# ============================================================
# CONFIG -- adjust these paths to match your directory layout
# ============================================================

BASE_DIR = Path(__file__).parent

IDL_MC_DIR = BASE_DIR / "data" / "all_IDL_MC_synthetic_dem_data"
IDL_MC_PATTERN = "idl_synthetic_dem_idx{n}_allfilters_MC1000.sav"

XRTPY_MC_DIR = BASE_DIR / "data" / "output_data" / "xrtpy_monte_carlo_dem_results"
XRTPY_MC_PATTERN = "xrtpy_dem_idx{n}_allfilters_MC1000.npz"

TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
TRUE_DEM_PATTERN = "DEM_{n}.txt"

CHI_SQUARE_TABLE_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "chi_square_comparison_tables"
)

OUT_LONG_CSV = (
    CHI_SQUARE_TABLE_DIR
    / "mc1000_chi_square_all_runs.csv"
)

OUT_SUMMARY_CSV = (
    CHI_SQUARE_TABLE_DIR
    / "mc1000_chi_square_summary.csv"
)

PLOT_DIR = (
    BASE_DIR
    / "plots"
    / "idl_xrtpy_synthetic_monte_carlo_1000_run_comparison"
)

OVERLAY_DIR = (
    PLOT_DIR
    / "per_dem_monte_carlo_overlays"
)

CHI_SQUARE_PLOT_DIR = (
    PLOT_DIR
    / "chi_square_summary_plots"
)

MEDIAN_HIST_PNG = (
    CHI_SQUARE_PLOT_DIR
    / "mc1000_median_chi_square_distribution_histogram.png"
)

ALLRUNS_HIST_PNG = (
    CHI_SQUARE_PLOT_DIR
    / "mc1000_all_runs_chi_square_distribution_histogram.png"
)

DEM_FLOOR = 1e-99
N_OVERLAY_MAX_CURVES = 1000  # set lower (e.g. 200) to speed up plotting if desired


# ============================================================
# Discovery (same pattern as compare_base_dem.py)
# ============================================================

def _extract_index(filename: str, pattern: str) -> int | None:
    regex = re.escape(pattern).replace(r"\{n\}", r"(\d+)")
    m = re.match(regex, filename)
    return int(m.group(1)) if m else None



def discover_common_indices() -> list[int]:
    """Find DEM indices available in both the IDL and XRTpy MC result sets."""
    idl_indices = {
        _extract_index(path.name, IDL_MC_PATTERN)
        for path in IDL_MC_DIR.glob(IDL_MC_PATTERN.format(n="*"))
    }
    idl_indices.discard(None)

    xrtpy_indices = {
        _extract_index(path.name, XRTPY_MC_PATTERN)
        for path in XRTPY_MC_DIR.glob(XRTPY_MC_PATTERN.format(n="*"))
    }
    xrtpy_indices.discard(None)

    common_indices = sorted(idl_indices & xrtpy_indices)
    idl_only = sorted(idl_indices - xrtpy_indices)
    xrtpy_only = sorted(xrtpy_indices - idl_indices)

    print(f"IDL MC files found:       {len(idl_indices)}")
    print(f"XRTpy MC files found:     {len(xrtpy_indices)}")
    print(f"Common indices (usable):  {len(common_indices)}")

    if idl_only:
        print(f"  IDL-only (missing XRTpy): {idl_only}")

    if xrtpy_only:
        print(f"  XRTpy-only (missing IDL): {xrtpy_only}")

    return common_indices


# ============================================================
# Loaders
# ============================================================


def load_idl_mc(dem_id: int) -> dict:
    """Load one IDL MC1000 result."""
    path = IDL_MC_DIR / IDL_MC_PATTERN.format(n=dem_id)

    data = readsav(str(path), python_dict=True)

    dem_out = np.asarray(data["dem_out"], dtype=float)
    chisq = np.asarray(data["chisq"], dtype=float).ravel()

    return {
        "logT": np.asarray(data["logt_out"], dtype=float).ravel(),
        "dem_base": dem_out[0, :],
        "dem_mc": dem_out[1:, :],
        "chisq_base": float(chisq[0]),
        "chisq_mc": chisq[1:],
    }


def load_xrtpy_mc(dem_id: int) -> dict:
    """Load one XRTpy MC1000 result."""
    path = XRTPY_MC_DIR / XRTPY_MC_PATTERN.format(n=dem_id)

    with np.load(path) as result:
        mc_dem = np.asarray(result["mc_dem"], dtype=float)
        mc_chisq = np.asarray(
            result["mc_chisq"],
            dtype=float,
        ).ravel()
        log_temperature = np.asarray(
            result["logT"],
            dtype=float,
        ).ravel()

    return {
        "logT": log_temperature,
        "dem_base": mc_dem[0, :],
        "dem_mc": mc_dem[1:, :],
        "chisq_base": float(mc_chisq[0]),
        "chisq_mc": mc_chisq[1:],
    }


def load_true_dem(dem_id: int) -> dict:
    """Load the synthetic truth DEM for one case."""
    path = TRUE_DEM_DIR / TRUE_DEM_PATTERN.format(n=dem_id)

    data = np.loadtxt(path)

    return {
        "logT": data[:, 0].astype(float),
        "dem": data[:, 1].astype(float),
    }


def _log10_dem(
    dem: np.ndarray,
    floor: float = DEM_FLOOR,
) -> np.ndarray:
    """Convert linear DEM values to log10 with a numerical floor."""
    return np.log10(
        np.maximum(
            np.asarray(dem, dtype=float),
            floor,
        )
    )


# ============================================================
# Step 2: MC ensemble overlay plots
# ============================================================
def plot_mc_overlay(
    dem_id: int,
    outdir: Path = OVERLAY_DIR,
) -> None:
    """Plot IDL and XRTpy MC ensembles against the synthetic truth."""
    outdir.mkdir(parents=True, exist_ok=True)

    idl = load_idl_mc(dem_id)
    xrtpy = load_xrtpy_mc(dem_id)
    truth = load_true_dem(dem_id)

    if not np.allclose(
        idl["logT"],
        xrtpy["logT"],
        atol=1e-8,
    ):
        print(
            f"  [warn] DEM {dem_id}: "
            "IDL and XRTpy logT grids differ."
        )

    log_temperature = xrtpy["logT"]

    fig, ax = plt.subplots(figsize=(9, 6.5))

    n_curves = min(
        N_OVERLAY_MAX_CURVES,
        idl["dem_mc"].shape[0],
        xrtpy["dem_mc"].shape[0],
    )

    for mc_index in range(n_curves):
        ax.step(
            log_temperature,
            _log10_dem(idl["dem_mc"][mc_index]),
            where="mid",
            color="orange",
            alpha=0.03,
            linewidth=0.8,
        )

    for mc_index in range(n_curves):
        ax.step(
            log_temperature,
            _log10_dem(xrtpy["dem_mc"][mc_index]),
            where="mid",
            color="steelblue",
            alpha=0.03,
            linewidth=0.8,
        )

    ax.step(
        truth["logT"],
        _log10_dem(truth["dem"]),
        where="mid",
        color="red",
        linewidth=2.4,
        label="Synthetic DEM (truth)",
    )

    ax.step(
        log_temperature,
        _log10_dem(idl["dem_base"]),
        where="mid",
        color="darkorange",
        linewidth=2.2,
        label="IDL base + MC",
    )

    ax.step(
        log_temperature,
        _log10_dem(xrtpy["dem_base"]),
        where="mid",
        color="navy",
        linestyle="--",
        linewidth=2.2,
        label="XRTpy base + MC",
    )

    idl_median_chisq = np.median(idl["chisq_mc"])
    xrtpy_median_chisq = np.median(xrtpy["chisq_mc"])

    ax.set_xlabel(r"$\log_{10} T$  [K]")
    ax.set_ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")

    ax.set_title(
        f"DEM {dem_id}  (MC1000)  |  "
        f"IDL median $\\chi^2$={idl_median_chisq:.3g}   "
        f"XRTpy median $\\chi^2$={xrtpy_median_chisq:.3g}"
    )

    ax.grid(alpha=0.3)
    ax.legend(loc="best")

    fig.tight_layout()

    out_png = outdir / f"DEM_{dem_id}_mc_overlay.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)



# ============================================================
# Step 4: Chi-square tables + summary plots
# ============================================================
def build_chisq_tables(
    indices: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Build per-run and per-DEM Monte Carlo chi-square comparison tables.

    Parameters
    ----------
    indices : list[int]
        DEM indices available in both the IDL and XRTpy MC result sets.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        The first table contains one row per paired Monte Carlo realization.
        The second contains per-DEM chi-square summary statistics.
    """
    long_rows = []
    summary_rows = []

    for dem_id in indices:
        try:
            idl = load_idl_mc(dem_id)
            xrtpy = load_xrtpy_mc(dem_id)
        except FileNotFoundError as error:
            print(f"  [skip] DEM {dem_id}: {error}")
            continue

        n_runs = min(
            len(idl["chisq_mc"]),
            len(xrtpy["chisq_mc"]),
        )

        for run_index in range(n_runs):
            long_rows.append(
                {
                    "dem_id": dem_id,
                    "run": run_index + 1,
                    "chisq_idl": idl["chisq_mc"][run_index],
                    "chisq_xrtpy": xrtpy["chisq_mc"][run_index],
                }
            )

        summary_rows.append(
            {
                "dem_id": dem_id,
                "chisq_idl_median": np.median(idl["chisq_mc"]),
                "chisq_idl_p16": np.percentile(
                    idl["chisq_mc"],
                    16,
                ),
                "chisq_idl_p84": np.percentile(
                    idl["chisq_mc"],
                    84,
                ),
                "chisq_xrtpy_median": np.median(
                    xrtpy["chisq_mc"]
                ),
                "chisq_xrtpy_p16": np.percentile(
                    xrtpy["chisq_mc"],
                    16,
                ),
                "chisq_xrtpy_p84": np.percentile(
                    xrtpy["chisq_mc"],
                    84,
                ),
            }
        )

    long_df = pd.DataFrame(long_rows)

    summary_df = (
        pd.DataFrame(summary_rows)
        .sort_values("dem_id")
        .reset_index(drop=True)
    )

    return long_df, summary_df



def plot_median_chisq_histogram(
    summary_df: pd.DataFrame,
    out_png: Path = MEDIAN_HIST_PNG,
) -> None:
    """Plot the distribution of per-DEM median Monte Carlo chi-square values."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    median_chisq_idl = summary_df[
        "chisq_idl_median"
    ].to_numpy(dtype=float)

    median_chisq_xrtpy = summary_df[
        "chisq_xrtpy_median"
    ].to_numpy(dtype=float)

    log_median_chisq_idl = np.log10(
        np.maximum(median_chisq_idl, 1e-10)
    )

    log_median_chisq_xrtpy = np.log10(
        np.maximum(median_chisq_xrtpy, 1e-10)
    )

    idl_median_of_medians = np.median(median_chisq_idl)
    xrtpy_median_of_medians = np.median(median_chisq_xrtpy)

    bins = np.linspace(
        min(
            log_median_chisq_idl.min(),
            log_median_chisq_xrtpy.min(),
        ),
        max(
            log_median_chisq_idl.max(),
            log_median_chisq_xrtpy.max(),
        ),
        30,
    )

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.hist(
        log_median_chisq_idl,
        bins=bins,
        alpha=0.6,
        color="orange",
        label="IDL",
    )

    ax.hist(
        log_median_chisq_xrtpy,
        bins=bins,
        alpha=0.6,
        color="steelblue",
        label="XRTpy",
    )

    ax.axvline(
        np.log10(max(idl_median_of_medians, 1e-10)),
        color="darkorange",
        linestyle="--",
        label=(
            "IDL median-of-medians "
            f"$\\chi^2$ = {idl_median_of_medians:.3g}"
        ),
    )

    ax.axvline(
        np.log10(max(xrtpy_median_of_medians, 1e-10)),
        color="navy",
        linestyle="--",
        label=(
            "XRTpy median-of-medians "
            f"$\\chi^2$ = {xrtpy_median_of_medians:.3g}"
        ),
    )

    ax.set_title(
        f"Median MC $\\chi^2$ Distribution — "
        f"{len(summary_df)} Synthetic DEMs "
        f"({N_OVERLAY_MAX_CURVES} MC each)"
    )

    ax.set_xlabel(
        r"$\log_{10}$(median $\chi^2$)  (per DEM)"
    )
    ax.set_ylabel("Number of DEMs")

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_all_runs_chisq_histogram(
    long_df: pd.DataFrame,
    out_png: Path = ALLRUNS_HIST_PNG,
) -> None:
    """Plot pooled chi-square distributions for all MC realizations."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    chisq_idl = long_df["chisq_idl"].to_numpy(dtype=float)
    chisq_xrtpy = long_df["chisq_xrtpy"].to_numpy(dtype=float)

    log_chisq_idl = np.log10(
        np.maximum(chisq_idl, 1e-10)
    )
    log_chisq_xrtpy = np.log10(
        np.maximum(chisq_xrtpy, 1e-10)
    )

    bins = np.linspace(
        min(
            log_chisq_idl.min(),
            log_chisq_xrtpy.min(),
        ),
        max(
            log_chisq_idl.max(),
            log_chisq_xrtpy.max(),
        ),
        60,
    )

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.hist(
        log_chisq_idl,
        bins=bins,
        alpha=0.5,
        color="orange",
        label=f"IDL (N={len(log_chisq_idl)} runs)",
    )

    ax.hist(
        log_chisq_xrtpy,
        bins=bins,
        alpha=0.5,
        color="steelblue",
        label=f"XRTpy (N={len(log_chisq_xrtpy)} runs)",
    )

    n_dems = long_df["dem_id"].nunique()

    ax.set_title(
        f"All Individual MC Run $\\chi^2$ — "
        f"{n_dems} Synthetic DEMs"
    )
    ax.set_xlabel(
        r"$\log_{10}(\chi^2)$  (single MC run)"
    )
    ax.set_ylabel("Number of runs")

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Main
def main():
    """Run the full Monte Carlo DEM comparison workflow."""
    print("Discovering usable DEM indices (MC1000)...")
    indices = discover_common_indices()

    if not indices:
        raise RuntimeError(
            "No common IDL and XRTpy MC1000 DEM results were found. "
            "Check IDL_MC_DIR and XRTPY_MC_DIR."
        )

    print("\nBuilding chi-square tables (long + summary)...")
    long_df, summary_df = build_chisq_tables(indices)

    if long_df.empty or summary_df.empty:
        raise RuntimeError(
            "No valid paired Monte Carlo comparisons could be built "
            "from the discovered result files."
        )

    OUT_LONG_CSV.parent.mkdir(parents=True, exist_ok=True)

    long_df.to_csv(
        OUT_LONG_CSV,
        index=False,
    )
    summary_df.to_csv(
        OUT_SUMMARY_CSV,
        index=False,
    )

    print(
        f"Saved: {OUT_LONG_CSV} "
        f"({len(long_df)} rows)"
    )
    print(
        f"Saved: {OUT_SUMMARY_CSV} "
        f"({len(summary_df)} rows)"
    )

    print("\nGenerating per-DEM MC overlay plots...")

    for dem_id in summary_df["dem_id"]:
        plot_mc_overlay(int(dem_id))

    print(
        f"Saved {len(summary_df)} overlay plots to: "
        f"{OVERLAY_DIR}"
    )

    print("\nGenerating chi-square summary plots...")

    plot_median_chisq_histogram(summary_df)
    plot_all_runs_chisq_histogram(long_df)

    print(f"Saved: {MEDIAN_HIST_PNG}")
    print(f"Saved: {ALLRUNS_HIST_PNG}")

    print("\nSummary:")
    print(f"  N DEMs = {len(summary_df)}")
    print(f"  N total MC runs compared = {len(long_df)}")
    print(
        "  IDL   median-of-medians chisq = "
        f"{summary_df['chisq_idl_median'].median():.4g}"
    )
    print(
        "  XRTpy median-of-medians chisq = "
        f"{summary_df['chisq_xrtpy_median'].median():.4g}"
    )


if __name__ == "__main__":
    main()

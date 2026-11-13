"""
Compare synthetic truth, IDL base DEM results, and XRTpy base DEM results.

The script performs two parts of the synthetic validation workflow:

1. Compare the recovered base DEM from IDL and XRTpy against the known
   synthetic truth for every DEM index available in both result sets.
2. Compare the corresponding base chi-square values across the campaign.

The script automatically discovers DEM indices that have both an IDL base SAV
file and an XRTpy base NPZ file, so partial result sets are handled without a
hardcoded index list.

Inputs
------
IDL base results:
    data/all_IDL_synthetic_dem_data/
    idl_synthetic_dem_idx{N}_allfilters_base.sav

XRTpy base results:
    data/output_data/xrtpy_base_dem_results/
    xrtpy_dem_idx{N}_allfilters_base.npz

Synthetic truth DEMs:
    data/synthetic_dems_data/
    DEM_{N}.txt

Outputs
-------
Comparison table:
    data/output_data/chi_square_comparison_tables/base_dem_comparison.csv

Per-DEM overlays:
    plots/idl_xrtpy_synthetic_base_dem_comparison/per_dem_overlays/

Chi-square summary plots:
    plots/idl_xrtpy_synthetic_base_dem_comparison/chi_square_summary_plots/
"""

from __future__ import annotations

import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.io import readsav


# ============================================================
# Configuration
# ============================================================

BASE_DIR = Path(__file__).parent

IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
IDL_BASE_PATTERN = "idl_synthetic_dem_idx{n}_allfilters_base.sav"

XRTPY_BASE_DIR = (
    BASE_DIR
    / "data"
    / "output_data"
    / "xrtpy_base_dem_results"
)
XRTPY_BASE_PATTERN = "xrtpy_dem_idx{n}_allfilters_base.npz"

TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
TRUE_DEM_PATTERN = "DEM_{n}.txt"

OUT_CSV = (
    BASE_DIR
    / "data"
    / "output_data"
    / "chi_square_comparison_tables"
    / "base_dem_comparison.csv"
)

PLOT_DIR = (
    BASE_DIR
    / "plots"
    / "idl_xrtpy_synthetic_base_dem_comparison"
)

OVERLAY_DIR = PLOT_DIR / "per_dem_overlays"
CHI_SQUARE_PLOT_DIR = PLOT_DIR / "chi_square_summary_plots"

HIST_PNG = (
    CHI_SQUARE_PLOT_DIR
    / "base_chi_square_distribution_histogram.png"
)

SCATTER_PNG = (
    CHI_SQUARE_PLOT_DIR
    / "base_chi_square_idl_xrtpy_scatter.png"
)

DEM_FLOOR = 1e-99

# ============================================================
# Discovery
# ============================================================

def _extract_index(filename: str, pattern: str) -> int | None:
    """Turn a glob-matched filename into its DEM index using `pattern`."""
    regex = re.escape(pattern).replace(r"\{n\}", r"(\d+)")
    m = re.match(regex, filename)
    return int(m.group(1)) if m else None


def discover_common_indices() -> list[int]:
    """
    Find DEM indices that have BOTH an IDL base .sav and an XRTpy base .npz.
    This is deliberately auto-discovering rather than hardcoded 0..161,
    since IDL runs may skip cases (nosolve, etc.).
    """
    idl_indices = set()
    for f in IDL_BASE_DIR.glob(IDL_BASE_PATTERN.format(n="*")):
        idx = _extract_index(f.name, IDL_BASE_PATTERN)
        if idx is not None:
            idl_indices.add(idx)

    xrtpy_indices = set()
    for f in XRTPY_BASE_DIR.glob(XRTPY_BASE_PATTERN.format(n="*")):
        idx = _extract_index(f.name, XRTPY_BASE_PATTERN)
        if idx is not None:
            xrtpy_indices.add(idx)

    common = sorted(idl_indices & xrtpy_indices)
    idl_only = sorted(idl_indices - xrtpy_indices)
    xrtpy_only = sorted(xrtpy_indices - idl_indices)

    print(f"IDL base files found:    {len(idl_indices)}")
    print(f"XRTpy base files found:  {len(xrtpy_indices)}")
    print(f"Common indices (usable): {len(common)}")
    if idl_only:
        print(f"  IDL-only (missing XRTpy):  {idl_only}")
    if xrtpy_only:
        print(f"  XRTpy-only (missing IDL):  {xrtpy_only}")

    return common


# ============================================================
# Loaders
# ============================================================

def load_idl_base(idx: int) -> dict:
    path = IDL_BASE_DIR / IDL_BASE_PATTERN.format(n=idx)
    d = readsav(str(path), python_dict=True)
    return {
        "logT": np.asarray(d["logt_out"], dtype=float).ravel(),
        "dem": np.asarray(d["dem_out"], dtype=float).ravel(),
        "chisq": float(np.asarray(d["chisq"]).ravel()[0]),
    }


def load_xrtpy_base(idx: int) -> dict:
    path = XRTPY_BASE_DIR / XRTPY_BASE_PATTERN.format(n=idx)
    d = np.load(path)
    return {
        "logT": np.asarray(d["logT"], dtype=float).ravel(),
        "dem": np.asarray(d["dem"], dtype=float).ravel(),
        "chisq": float(np.asarray(d["chisq"]).ravel()[0]),
    }


def load_true_dem(idx: int) -> dict:
    path = TRUE_DEM_DIR / TRUE_DEM_PATTERN.format(n=idx)
    arr = np.loadtxt(path)
    return {
        "logT": arr[:, 0].astype(float),
        "dem": arr[:, 1].astype(float),
    }


def _log10_dem(dem: np.ndarray, floor: float = DEM_FLOOR) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))


# ============================================================
# Step 1: Base DEM comparison (build table + per-DEM plots)
# ============================================================


def build_comparison_table(indices: list[int]) -> pd.DataFrame:
    """Build the paired IDL/XRTpy base chi-square comparison table."""
    rows = []

    for dem_id in indices:
        try:
            idl = load_idl_base(dem_id)
            xrtpy = load_xrtpy_base(dem_id)

            # Load the truth DEM to ensure the corresponding synthetic case
            # is present before including it in the comparison.
            load_true_dem(dem_id)

        except FileNotFoundError as error:
            print(f"  [skip] DEM {dem_id}: {error}")
            continue

        if not np.allclose(
            idl["logT"],
            xrtpy["logT"],
            atol=1e-8,
        ):
            print(
                f"  [warn] DEM {dem_id}: "
                "IDL and XRTpy logT grids differ."
            )

        rows.append(
            {
                "dem_id": dem_id,
                "chisq_idl": idl["chisq"],
                "chisq_xrtpy": xrtpy["chisq"],
                "delta_chisq": xrtpy["chisq"] - idl["chisq"],
            }
        )

    return (
        pd.DataFrame(rows)
        .sort_values("dem_id")
        .reset_index(drop=True)
    )



def plot_base_overlay(
    dem_id: int,
    outdir: Path = OVERLAY_DIR,
) -> None:
    """Plot synthetic truth, IDL base DEM, and XRTpy base DEM."""
    outdir.mkdir(parents=True, exist_ok=True)

    truth = load_true_dem(dem_id)
    idl = load_idl_base(dem_id)
    xrtpy = load_xrtpy_base(dem_id)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.step(
        truth["logT"],
        _log10_dem(truth["dem"]),
        where="mid",
        color="red",
        linewidth=2.2,
        label="Synthetic DEM (truth)",
    )

    ax.step(
        idl["logT"],
        _log10_dem(idl["dem"]),
        where="mid",
        color="orange",
        linewidth=2.0,
        label="IDL base",
    )

    ax.step(
        xrtpy["logT"],
        _log10_dem(xrtpy["dem"]),
        where="mid",
        color="blue",
        linestyle="--",
        linewidth=2.0,
        label="XRTpy base",
    )

    ax.set_xlabel(r"$\log_{10} T$  [K]")
    ax.set_ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")

    ax.set_title(
        f"DEM {dem_id}  |  "
        f"IDL $\\chi^2$={idl['chisq']:.3g}   "
        f"XRTpy $\\chi^2$={xrtpy['chisq']:.3g}"
    )

    ax.grid(alpha=0.3)
    ax.legend()

    fig.tight_layout()

    out_png = outdir / f"DEM_{dem_id}_base_overlay.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Step 3: Chi-square summary (histogram + paired scatter)
# ============================================================

def plot_chisq_histogram(
    df: pd.DataFrame,
    out_png: Path = HIST_PNG,
) -> None:
    """Plot the IDL and XRTpy base chi-square distributions."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    chisq_idl = df["chisq_idl"].to_numpy(dtype=float)
    chisq_xrtpy = df["chisq_xrtpy"].to_numpy(dtype=float)

    log_chisq_idl = np.log10(np.maximum(chisq_idl, 1e-10))
    log_chisq_xrtpy = np.log10(np.maximum(chisq_xrtpy, 1e-10))

    median_idl = np.median(chisq_idl)
    median_xrtpy = np.median(chisq_xrtpy)

    bins = np.linspace(
        min(log_chisq_idl.min(), log_chisq_xrtpy.min()),
        max(log_chisq_idl.max(), log_chisq_xrtpy.max()),
        30,
    )

    fig, ax = plt.subplots(figsize=(9, 6))

    ax.hist(
        log_chisq_idl,
        bins=bins,
        alpha=0.6,
        color="orange",
        label="IDL",
    )
    ax.hist(
        log_chisq_xrtpy,
        bins=bins,
        alpha=0.6,
        color="steelblue",
        label="XRTpy",
    )

    ax.axvline(
        np.log10(max(median_idl, 1e-10)),
        color="darkorange",
        linestyle="--",
        label=f"IDL median $\\chi^2$ = {median_idl:.3g}",
    )
    ax.axvline(
        np.log10(max(median_xrtpy, 1e-10)),
        color="navy",
        linestyle="--",
        label=f"XRTpy median $\\chi^2$ = {median_xrtpy:.3g}",
    )

    ax.set_title(
        f"Base DEM $\\chi^2$ Distribution — "
        f"{len(df)} Synthetic DEMs"
    )
    ax.set_xlabel(r"$\log_{10}(\chi^2)$  (base DEM)")
    ax.set_ylabel("Number of DEMs")

    ax.legend()
    ax.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)



def plot_chisq_scatter(
    df: pd.DataFrame,
    out_png: Path = SCATTER_PNG,
) -> None:
    """Plot paired IDL and XRTpy base chi-square values by IDL rank."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    ranked = (
        df.sort_values("chisq_idl")
        .reset_index(drop=True)
    )
    rank = np.arange(len(ranked))

    fig, ax = plt.subplots(figsize=(11, 6))

    ax.scatter(
        rank,
        ranked["chisq_idl"],
        color="orange",
        label="IDL",
        s=25,
        alpha=0.8,
    )

    ax.scatter(
        rank,
        ranked["chisq_xrtpy"],
        color="steelblue",
        label="XRTpy",
        s=25,
        alpha=0.8,
    )

    ax.set_yscale("log")
    ax.set_xlabel(r"DEM rank (sorted by IDL $\chi^2$)")
    ax.set_ylabel(r"$\chi^2$ (base DEM, log scale)")
    ax.set_title(
        f"Base DEM $\\chi^2$ per Case — "
        f"{len(ranked)} Synthetic DEMs"
    )

    ax.legend()
    ax.grid(alpha=0.3, which="both")

    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Main
# ============================================================
def main():
    """Run the full base DEM comparison workflow."""
    print("Discovering usable DEM indices...")
    indices = discover_common_indices()

    if not indices:
        raise RuntimeError(
            "No common IDL and XRTpy base DEM results were found. "
            "Check IDL_BASE_DIR and XRTPY_BASE_DIR."
        )

    print("\nBuilding chi-square comparison table...")
    comparison = build_comparison_table(indices)

    if comparison.empty:
        raise RuntimeError(
            "No valid base DEM comparisons could be built from the "
            "discovered result files."
        )

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    comparison.to_csv(OUT_CSV, index=False)

    print(f"Saved: {OUT_CSV} ({len(comparison)} rows)")

    print("\nGenerating per-DEM overlay plots...")

    for dem_id in comparison["dem_id"]:
        plot_base_overlay(int(dem_id))

    print(
        f"Saved {len(comparison)} overlay plots to: "
        f"{OVERLAY_DIR}"
    )

    print("\nGenerating chi-square summary plots...")

    plot_chisq_histogram(comparison)
    plot_chisq_scatter(comparison)

    print(f"Saved: {HIST_PNG}")
    print(f"Saved: {SCATTER_PNG}")

    n_xrtpy_lower = (
        comparison["delta_chisq"] < 0
    ).sum()

    print("\nSummary:")
    print(f"  N = {len(comparison)}")
    print(
        "  IDL   median chisq = "
        f"{comparison['chisq_idl'].median():.4g}"
    )
    print(
        "  XRTpy median chisq = "
        f"{comparison['chisq_xrtpy'].median():.4g}"
    )
    print(
        "  Cases where XRTpy chisq < IDL chisq: "
        f"{n_xrtpy_lower} / {len(comparison)}"
    )


if __name__ == "__main__":
    main()
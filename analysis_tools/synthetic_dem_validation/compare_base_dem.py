"""
compare_base_dem.py

Base DEM comparison pipeline (#1 and #3 from the validation plan):

  1. Base DEM comparison: Synthetic (truth) vs. IDL vs. XRTpy
  3. Chi-square of every base DEM, IDL vs. XRTpy

Auto-discovers which DEM indices have BOTH an IDL base .sav and an
XRTpy base .npz, so partial runs (e.g. IDL's 139/162 solved) are
handled without manual index lists.

Expected inputs (adjust CONFIG block below if your paths differ):

  IDL base .sav   : data/all_IDL_synthetic_dem_data/idl_synthetic_dem_idx{N}_allfilters_base.sav
      keys used: logt_out, dem_out, chisq

  XRTpy base .npz : xrtpy_output/xrtpy_dem_idx{N}_allfilters_base.npz
      keys used: logT, dem, chisq

  True synthetic DEM (ground truth), two-column logT/linear DEM:
      synthetic_dems_data/DEM_{N}.txt
      (assumed whitespace-delimited: logT  DEM)

Outputs:
  data/base_dem_comparison.csv   -- one row per DEM index with both chisq values
  plots/base_dem_overlay/DEM_{N}_base_overlay.png  -- per-DEM 3-curve plot
  plots/base_chisq_histogram.png                    -- summary histogram
  plots/base_chisq_scatter.png                      -- paired-marker plot (IDL vs XRTpy)
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

BASE_DIR = Path(".")

IDL_BASE_DIR = BASE_DIR / "data" / "all_IDL_synthetic_dem_data"
IDL_BASE_PATTERN = "idl_synthetic_dem_idx{n}_allfilters_base.sav"

XRTPY_BASE_DIR = BASE_DIR / "xrtpy_output"
XRTPY_BASE_PATTERN = "xrtpy_dem_idx{n}_allfilters_base.npz"

TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
TRUE_DEM_PATTERN = "DEM_{n}.txt"

OUT_CSV = BASE_DIR / "data" / "base_dem_comparison.csv"
OVERLAY_DIR = BASE_DIR / "plots" / "base_dem_overlay"
HIST_PNG = BASE_DIR / "plots" / "base_chisq_histogram.png"
SCATTER_PNG = BASE_DIR / "plots" / "base_chisq_scatter.png"

DEM_FLOOR = 1e-99  # for log10 safety


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
    rows = []
    for idx in indices:
        try:
            idl = load_idl_base(idx)
            xrt = load_xrtpy_base(idx)
            truth = load_true_dem(idx)
        except FileNotFoundError as e:
            print(f"  [skip] DEM {idx}: {e}")
            continue

        if not np.allclose(idl["logT"], xrt["logT"], atol=1e-8):
            print(f"  [warn] DEM {idx}: IDL and XRTpy logT grids differ!")

        rows.append(
            {
                "dem_id": idx,
                "chisq_idl": idl["chisq"],
                "chisq_xrtpy": xrt["chisq"],
                "delta_chisq": xrt["chisq"] - idl["chisq"],
            }
        )

    df = pd.DataFrame(rows).sort_values("dem_id").reset_index(drop=True)
    return df


def plot_base_overlay(idx: int, outdir: Path = OVERLAY_DIR) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    idl = load_idl_base(idx)
    xrt = load_xrtpy_base(idx)
    truth = load_true_dem(idx)

    fig, ax = plt.subplots(figsize=(8, 5.5))

    ax.step(
        truth["logT"], _log10_dem(truth["dem"]),
        where="mid", color="red", linewidth=2.2, label="Synthetic DEM (truth)",
    )
    ax.step(
        idl["logT"], _log10_dem(idl["dem"]),
        where="mid", color="orange", linewidth=2.0, label="IDL base",
    )
    ax.step(
        xrt["logT"], _log10_dem(xrt["dem"]),
        where="mid", color="blue", linestyle="--", linewidth=2.0, label="XRTpy base",
    )

    ax.set_xlabel(r"$\log_{10} T$  [K]")
    ax.set_ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    ax.set_title(
        f"DEM {idx}  |  IDL $\\chi^2$={idl['chisq']:.3g}   "
        f"XRTpy $\\chi^2$={xrt['chisq']:.3g}"
    )
    ax.grid(alpha=0.3)
    ax.legend()
    fig.tight_layout()

    out_png = outdir / f"DEM_{idx}_base_overlay.png"
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Step 3: Chi-square summary (histogram + paired scatter)
# ============================================================

def plot_chisq_histogram(df: pd.DataFrame, out_png: Path = HIST_PNG) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    log_idl = np.log10(np.maximum(df["chisq_idl"].to_numpy(), 1e-10))
    log_xrt = np.log10(np.maximum(df["chisq_xrtpy"].to_numpy(), 1e-10))

    fig, ax = plt.subplots(figsize=(9, 6))
    bins = np.linspace(
        min(log_idl.min(), log_xrt.min()),
        max(log_idl.max(), log_xrt.max()),
        30,
    )
    ax.hist(log_idl, bins=bins, alpha=0.6, color="orange", label="IDL")
    ax.hist(log_xrt, bins=bins, alpha=0.6, color="steelblue", label="XRTpy")
    ax.axvline(np.median(log_idl), color="darkorange", linestyle="--",
               label=f"IDL median $\\chi^2$ = {10**np.median(log_idl):.3g}")
    ax.axvline(np.median(log_xrt), color="navy", linestyle="--",
               label=f"XRTpy median $\\chi^2$ = {10**np.median(log_xrt):.3g}")

    ax.set_title(f"Base DEM $\\chi^2$ Distribution — {len(df)} synthetic DEMs")
    ax.set_xlabel(r"$\log_{10}(\chi^2)$  (base DEM)")
    ax.set_ylabel("Number of DEMs")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_chisq_scatter(df: pd.DataFrame, out_png: Path = SCATTER_PNG) -> None:
    """Paired-marker plot: DEMs ranked by IDL chisq, IDL vs XRTpy shown side by side."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    d = df.sort_values("chisq_idl").reset_index(drop=True)
    x = np.arange(len(d))

    fig, ax = plt.subplots(figsize=(11, 6))
    ax.scatter(x, d["chisq_idl"], color="orange", label="IDL", s=25, alpha=0.8)
    ax.scatter(x, d["chisq_xrtpy"], color="steelblue", label="XRTpy", s=25, alpha=0.8)
    ax.set_yscale("log")
    ax.set_xlabel("DEM rank (sorted by IDL $\\chi^2$)")
    ax.set_ylabel(r"$\chi^2$ (base DEM, log scale)")
    ax.set_title(f"Base DEM $\\chi^2$ per case — {len(d)} synthetic DEMs")
    ax.legend()
    ax.grid(alpha=0.3, which="both")
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    print("Discovering usable DEM indices...")
    indices = discover_common_indices()

    print("\nBuilding chi-square comparison table...")
    df = build_comparison_table(indices)

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(OUT_CSV, index=False)
    print(f"Saved: {OUT_CSV}  ({len(df)} rows)")

    print("\nGenerating per-DEM overlay plots...")
    for idx in df["dem_id"]:
        plot_base_overlay(int(idx))
    print(f"Saved {len(df)} overlay plots to: {OVERLAY_DIR}")

    print("\nGenerating chi-square summary plots...")
    plot_chisq_histogram(df)
    plot_chisq_scatter(df)
    print(f"Saved: {HIST_PNG}")
    print(f"Saved: {SCATTER_PNG}")

    print("\nSummary:")
    print(f"  N = {len(df)}")
    print(f"  IDL   median chisq = {df['chisq_idl'].median():.4g}")
    print(f"  XRTpy median chisq = {df['chisq_xrtpy'].median():.4g}")
    print(f"  Cases where XRTpy chisq < IDL chisq: "
          f"{(df['delta_chisq'] < 0).sum()} / {len(df)}")


if __name__ == "__main__":
    main()
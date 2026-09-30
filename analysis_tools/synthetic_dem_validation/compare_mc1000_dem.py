"""
compare_mc1000_dem.py

Monte Carlo (#2) and per-run chi-square (#4) comparison:

  2. MC1000 DEM comparison: IDL vs. XRTpy Monte Carlo ensembles
  4. Chi-square of every one of the 1000 MC realizations, IDL vs. XRTpy

Both sides share the same row convention:
  row 0        = base (unperturbed) solution
  rows 1..1000 = Monte Carlo perturbed realizations

  IDL   .sav : dem_out (1001, nT),  chisq (1001,)
  XRTpy .npz : mc_dem  (1001, nT),  mc_chisq (1001,)

Auto-discovers indices with BOTH an IDL MC .sav and an XRTpy MC .npz.

Outputs:
  data/mc1000_chisq_long.csv               -- one row per (dem_id, run) with both chisq
  data/mc1000_chisq_summary.csv            -- one row per dem_id: median/16/84 chisq for both solvers
  plots/mc1000_overlay/DEM_{N}_mc_overlay.png   -- per-DEM ensemble overlay (truth + IDL cloud + XRTpy cloud)
  plots/mc1000_chisq_histogram.png              -- median-MC-chisq histogram (per-DEM median, then pooled)
  plots/mc1000_chisq_all_runs_histogram.png     -- every single MC chisq pooled across all DEMs (not just medians)
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

IDL_MC_DIR = BASE_DIR / "data" / "all_IDL_MC_synthetic_dem_data"
IDL_MC_PATTERN = "idl_synthetic_dem_idx{n}_allfilters_MC1000.sav"

XRTPY_MC_DIR = BASE_DIR / "xrtpy_output"
XRTPY_MC_PATTERN = "xrtpy_dem_idx{n}_allfilters_MC1000.npz"

TRUE_DEM_DIR = BASE_DIR / "data" / "synthetic_dems_data"
TRUE_DEM_PATTERN = "DEM_{n}.txt"

OUT_LONG_CSV = BASE_DIR / "data" / "mc1000_chisq_long.csv"
OUT_SUMMARY_CSV = BASE_DIR / "data" / "mc1000_chisq_summary.csv"
OVERLAY_DIR = BASE_DIR / "plots" / "mc1000_overlay"
MEDIAN_HIST_PNG = BASE_DIR / "plots" / "mc1000_chisq_histogram.png"
ALLRUNS_HIST_PNG = BASE_DIR / "plots" / "mc1000_chisq_all_runs_histogram.png"

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
    idl_indices = {
        _extract_index(f.name, IDL_MC_PATTERN)
        for f in IDL_MC_DIR.glob(IDL_MC_PATTERN.format(n="*"))
    }
    idl_indices.discard(None)

    xrtpy_indices = {
        _extract_index(f.name, XRTPY_MC_PATTERN)
        for f in XRTPY_MC_DIR.glob(XRTPY_MC_PATTERN.format(n="*"))
    }
    xrtpy_indices.discard(None)

    common = sorted(idl_indices & xrtpy_indices)
    idl_only = sorted(idl_indices - xrtpy_indices)
    xrtpy_only = sorted(xrtpy_indices - idl_indices)

    print(f"IDL MC files found:    {len(idl_indices)}")
    print(f"XRTpy MC files found:  {len(xrtpy_indices)}")
    print(f"Common indices (usable): {len(common)}")
    if idl_only:
        print(f"  IDL-only (missing XRTpy):  {idl_only}")
    if xrtpy_only:
        print(f"  XRTpy-only (missing IDL):  {xrtpy_only}")

    return common


# ============================================================
# Loaders
# ============================================================

def load_idl_mc(idx: int) -> dict:
    path = IDL_MC_DIR / IDL_MC_PATTERN.format(n=idx)
    d = readsav(str(path), python_dict=True)
    dem_out = np.asarray(d["dem_out"], dtype=float)   # (1001, nT)
    chisq = np.asarray(d["chisq"], dtype=float).ravel()  # (1001,)
    return {
        "logT": np.asarray(d["logt_out"], dtype=float).ravel(),
        "dem_base": dem_out[0, :],
        "dem_mc": dem_out[1:, :],       # (1000, nT)
        "chisq_base": float(chisq[0]),
        "chisq_mc": chisq[1:],          # (1000,)
    }


def load_xrtpy_mc(idx: int) -> dict:
    path = XRTPY_MC_DIR / XRTPY_MC_PATTERN.format(n=idx)
    d = np.load(path)
    mc_dem = np.asarray(d["mc_dem"], dtype=float)      # (1001, nT)
    mc_chisq = np.asarray(d["mc_chisq"], dtype=float).ravel()  # (1001,)
    return {
        "logT": np.asarray(d["logT"], dtype=float).ravel(),
        "dem_base": mc_dem[0, :],
        "dem_mc": mc_dem[1:, :],
        "chisq_base": float(mc_chisq[0]),
        "chisq_mc": mc_chisq[1:],
    }


def load_true_dem(idx: int) -> dict:
    path = TRUE_DEM_DIR / TRUE_DEM_PATTERN.format(n=idx)
    arr = np.loadtxt(path)
    return {"logT": arr[:, 0].astype(float), "dem": arr[:, 1].astype(float)}


def _log10_dem(dem: np.ndarray, floor: float = DEM_FLOOR) -> np.ndarray:
    return np.log10(np.maximum(dem, floor))


# ============================================================
# Step 2: MC ensemble overlay plots
# ============================================================

def plot_mc_overlay(idx: int, outdir: Path = OVERLAY_DIR) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    idl = load_idl_mc(idx)
    xrt = load_xrtpy_mc(idx)
    truth = load_true_dem(idx)

    if not np.allclose(idl["logT"], xrt["logT"], atol=1e-8):
        print(f"  [warn] DEM {idx}: IDL and XRTpy logT grids differ!")

    logT = xrt["logT"]

    fig, ax = plt.subplots(figsize=(9, 6.5))

    n_curves = min(N_OVERLAY_MAX_CURVES, idl["dem_mc"].shape[0])

    # MC clouds (faint)
    for i in range(n_curves):
        ax.step(logT, _log10_dem(idl["dem_mc"][i]), where="mid",
                color="orange", alpha=0.03, linewidth=0.8)
    for i in range(n_curves):
        ax.step(logT, _log10_dem(xrt["dem_mc"][i]), where="mid",
                color="steelblue", alpha=0.03, linewidth=0.8)

    # Base curves (bold)
    ax.step(truth["logT"], _log10_dem(truth["dem"]), where="mid",
            color="red", linewidth=2.4, label="Synthetic DEM (truth)")
    ax.step(logT, _log10_dem(idl["dem_base"]), where="mid",
            color="darkorange", linewidth=2.2, label="IDL base + MC")
    ax.step(logT, _log10_dem(xrt["dem_base"]), where="mid",
            color="navy", linestyle="--", linewidth=2.2, label="XRTpy base + MC")

    ax.set_xlabel(r"$\log_{10} T$  [K]")
    ax.set_ylabel(r"$\log_{10}$ DEM  [cm$^{-5}$ K$^{-1}$]")
    ax.set_title(
        f"DEM {idx}  (MC1000)  |  IDL median $\\chi^2$={np.median(idl['chisq_mc']):.3g}   "
        f"XRTpy median $\\chi^2$={np.median(xrt['chisq_mc']):.3g}"
    )
    ax.grid(alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()

    out_png = outdir / f"DEM_{idx}_mc_overlay.png"
    fig.savefig(out_png, dpi=180)
    plt.close(fig)


# ============================================================
# Step 4: Chi-square tables + summary plots
# ============================================================

def build_chisq_tables(indices: list[int]) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Returns:
      long_df    -- one row per (dem_id, run 1..1000): chisq_idl, chisq_xrtpy
      summary_df -- one row per dem_id: median/16th/84th pct chisq for both solvers
    """
    long_rows = []
    summary_rows = []

    for idx in indices:
        try:
            idl = load_idl_mc(idx)
            xrt = load_xrtpy_mc(idx)
        except FileNotFoundError as e:
            print(f"  [skip] DEM {idx}: {e}")
            continue

        n_runs = min(len(idl["chisq_mc"]), len(xrt["chisq_mc"]))
        for run in range(n_runs):
            long_rows.append(
                {
                    "dem_id": idx,
                    "run": run + 1,
                    "chisq_idl": idl["chisq_mc"][run],
                    "chisq_xrtpy": xrt["chisq_mc"][run],
                }
            )

        summary_rows.append(
            {
                "dem_id": idx,
                "chisq_idl_median": np.median(idl["chisq_mc"]),
                "chisq_idl_p16": np.percentile(idl["chisq_mc"], 16),
                "chisq_idl_p84": np.percentile(idl["chisq_mc"], 84),
                "chisq_xrtpy_median": np.median(xrt["chisq_mc"]),
                "chisq_xrtpy_p16": np.percentile(xrt["chisq_mc"], 16),
                "chisq_xrtpy_p84": np.percentile(xrt["chisq_mc"], 84),
            }
        )

    long_df = pd.DataFrame(long_rows)
    summary_df = pd.DataFrame(summary_rows).sort_values("dem_id").reset_index(drop=True)
    return long_df, summary_df


def plot_median_chisq_histogram(summary_df: pd.DataFrame, out_png: Path = MEDIAN_HIST_PNG) -> None:
    out_png.parent.mkdir(parents=True, exist_ok=True)

    log_idl = np.log10(np.maximum(summary_df["chisq_idl_median"].to_numpy(), 1e-10))
    log_xrt = np.log10(np.maximum(summary_df["chisq_xrtpy_median"].to_numpy(), 1e-10))

    fig, ax = plt.subplots(figsize=(9, 6))
    bins = np.linspace(min(log_idl.min(), log_xrt.min()),
                        max(log_idl.max(), log_xrt.max()), 30)
    ax.hist(log_idl, bins=bins, alpha=0.6, color="orange", label="IDL")
    ax.hist(log_xrt, bins=bins, alpha=0.6, color="steelblue", label="XRTpy")
    ax.axvline(np.median(log_idl), color="darkorange", linestyle="--",
               label=f"IDL median $\\chi^2$ = {10**np.median(log_idl):.3g}")
    ax.axvline(np.median(log_xrt), color="navy", linestyle="--",
               label=f"XRTpy median $\\chi^2$ = {10**np.median(log_xrt):.3g}")

    ax.set_title(f"Median MC $\\chi^2$ Distribution — {len(summary_df)} synthetic DEMs (1000 MC each)")
    ax.set_xlabel(r"$\log_{10}$(median $\chi^2$)  (1000 MC per DEM)")
    ax.set_ylabel("Number of DEMs")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


def plot_all_runs_chisq_histogram(long_df: pd.DataFrame, out_png: Path = ALLRUNS_HIST_PNG) -> None:
    """Every single MC run's chisq, pooled across all DEMs (not just the per-DEM median)."""
    out_png.parent.mkdir(parents=True, exist_ok=True)

    log_idl = np.log10(np.maximum(long_df["chisq_idl"].to_numpy(), 1e-10))
    log_xrt = np.log10(np.maximum(long_df["chisq_xrtpy"].to_numpy(), 1e-10))

    fig, ax = plt.subplots(figsize=(9, 6))
    bins = np.linspace(min(log_idl.min(), log_xrt.min()),
                        max(log_idl.max(), log_xrt.max()), 60)
    ax.hist(log_idl, bins=bins, alpha=0.5, color="orange", label=f"IDL (N={len(log_idl)} runs)")
    ax.hist(log_xrt, bins=bins, alpha=0.5, color="steelblue", label=f"XRTpy (N={len(log_xrt)} runs)")

    ax.set_title(f"All Individual MC Run $\\chi^2$ — {long_df['dem_id'].nunique()} DEMs × 1000 MC")
    ax.set_xlabel(r"$\log_{10}(\chi^2)$  (single MC run)")
    ax.set_ylabel("Number of runs")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out_png, dpi=200)
    plt.close(fig)


# ============================================================
# Main
# ============================================================

def main():
    print("Discovering usable DEM indices (MC1000)...")
    indices = discover_common_indices()

    print("\nBuilding chi-square tables (long + summary)...")
    long_df, summary_df = build_chisq_tables(indices)

    OUT_LONG_CSV.parent.mkdir(parents=True, exist_ok=True)
    long_df.to_csv(OUT_LONG_CSV, index=False)
    summary_df.to_csv(OUT_SUMMARY_CSV, index=False)
    print(f"Saved: {OUT_LONG_CSV}  ({len(long_df)} rows)")
    print(f"Saved: {OUT_SUMMARY_CSV}  ({len(summary_df)} rows)")

    print("\nGenerating per-DEM MC overlay plots...")
    for idx in summary_df["dem_id"]:
        plot_mc_overlay(int(idx))
    print(f"Saved {len(summary_df)} overlay plots to: {OVERLAY_DIR}")

    print("\nGenerating chi-square summary plots...")
    plot_median_chisq_histogram(summary_df)
    plot_all_runs_chisq_histogram(long_df)
    print(f"Saved: {MEDIAN_HIST_PNG}")
    print(f"Saved: {ALLRUNS_HIST_PNG}")

    print("\nSummary:")
    print(f"  N DEMs = {len(summary_df)}")
    print(f"  N total MC runs compared = {len(long_df)}")
    print(f"  IDL   median-of-medians chisq = {summary_df['chisq_idl_median'].median():.4g}")
    print(f"  XRTpy median-of-medians chisq = {summary_df['chisq_xrtpy_median'].median():.4g}")


if __name__ == "__main__":
    main()
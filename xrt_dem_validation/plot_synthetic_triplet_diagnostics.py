
# """Additional visual diagnostics for 162 synthetic DEM inversions.

# Produces simple recovery plots, Python-IDL agreement diagnostics,
# and actual DEM curves. No inversions are performed.
# """

# from pathlib import Path

# import matplotlib.pyplot as plt
# import numpy as np
# import pandas as pd
# from scipy.io import readsav

# BASE = Path(__file__).resolve().parent
# DATA = BASE / "data"
# RESULTS = DATA / "output_data" / "scientific_validation"
# PLOTS = BASE / "plots" / "scientific_validation" / "diagnostics"

# PAIRS = [
#     ("python_vs_truth", "Python vs synthetic", "#2878b5"),
#     ("idl_vs_truth", "IDL vs synthetic", "#e68628"),
#     ("python_vs_idl", "Python vs IDL", "#666666"),
# ]


# def load_metrics():
#     df = pd.read_csv(RESULTS / "synthetic_triplet_metrics.csv")
#     df = df.sort_values("dem_id").reset_index(drop=True)

#     if len(df) != 162 or not np.array_equal(
#         df.dem_id.to_numpy(), np.arange(162)
#     ):
#         raise ValueError("Expected DEM IDs 0 through 161.")

#     return df


# def distance(df, pair, region="gray_band"):
#     return df[f"{pair}_{region}_distance"].to_numpy()


# def finish(fig, name):
#     fig.savefig(PLOTS / name, dpi=220, bbox_inches="tight")
#     plt.close(fig)
#     print("Saved:", name)


# def plot_simple_histograms(df):
#     """Easy-to-read distributions of the three comparisons."""
#     fig, ax = plt.subplots(figsize=(10, 6))

#     bins = np.linspace(0, 2, 21)

#     for pair, label, color in PAIRS:
#         values = distance(df, pair)
#         ax.hist(
#             values, bins=bins, histtype="step",
#             linewidth=2.2, color=color, label=label
#         )
#         ax.axvline(
#             np.median(values), color=color,
#             linestyle="--", linewidth=1.5
#         )

#     ax.set(
#         xlabel="Gray-band DEM discrepancy (0 = identical, 2 = maximum)",
#         ylabel="Number of synthetic DEM cases",
#         title="Recovery and implementation agreement: 162 DEMs",
#         xlim=(0, 2),
#     )
#     ax.legend(title="Dashed lines = medians")
#     ax.grid(alpha=0.2)

#     finish(fig, "01_simple_discrepancy_histograms.png")


# def plot_python_vs_idl_recovery(df):
#     """Directly compare recovery of synthetic truth case by case."""
#     py = distance(df, "python_vs_truth")
#     idl = distance(df, "idl_vs_truth")
#     delta = py - idl

#     fig, ax = plt.subplots(figsize=(8, 8))

#     colors = np.where(delta < -1e-10, "#2878b5", "#e68628")

#     ax.scatter(idl, py, c=colors, s=32, alpha=0.75)
#     ax.plot([0, 2], [0, 2], "k--", linewidth=1.3)

#     ax.text(
#         0.05, 0.95,
#         f"Python closer: {(delta < -1e-10).sum()}\n"
#         f"IDL closer: {(delta > 1e-10).sum()}\n"
#         f"Ties: {(np.abs(delta) <= 1e-10).sum()}",
#         transform=ax.transAxes, va="top",
#         bbox=dict(facecolor="white", alpha=0.9)
#     )

#     ax.set(
#         xlabel="IDL vs synthetic truth discrepancy",
#         ylabel="Python vs synthetic truth discrepancy",
#         title="Which solver reconstructs synthetic truth more closely?",
#         xlim=(0, 2.05),
#         ylim=(0, 2.05),
#     )
#     ax.set_aspect("equal")
#     ax.grid(alpha=0.2)

#     finish(fig, "02_python_vs_idl_truth_recovery.png")


# def plot_case_trends(df):
#     """Show discrepancies by DEM index to identify possible clusters."""
#     ids = df.dem_id.to_numpy()

#     fig, axes = plt.subplots(
#         3, 1, figsize=(13, 9), sharex=True, sharey=True
#     )

#     for ax, (pair, label, color) in zip(axes, PAIRS):
#         values = distance(df, pair)
#         ax.scatter(ids, values, s=22, color=color, alpha=0.75)
#         ax.plot(ids, values, color=color, alpha=0.2)

#         ax.set_ylabel("DEM discrepancy")
#         ax.set_title(label, loc="left")
#         ax.set_ylim(-0.05, 2.05)
#         ax.grid(alpha=0.2)

#     axes[-1].set_xlabel("Synthetic DEM ID")
#     fig.suptitle(
#         "Gray-band discrepancies across all synthetic DEM cases",
#         fontsize=14
#     )
#     fig.tight_layout()

#     finish(fig, "03_discrepancy_by_dem_id.png")


# def plot_agreement_heatmap(df):
#     """Compact overview of discrepancies for all cases."""
#     matrix = np.vstack([
#         distance(df, pair) for pair, _, _ in PAIRS
#     ])

#     fig, ax = plt.subplots(figsize=(15, 3.5))
#     im = ax.imshow(
#         matrix, aspect="auto", cmap="viridis",
#         vmin=0, vmax=2, interpolation="nearest"
#     )

#     ax.set_yticks(range(3))
#     ax.set_yticklabels([label for _, label, _ in PAIRS])
#     ax.set_xlabel("Synthetic DEM ID")
#     ax.set_title("Gray-band discrepancy map: lower is better")

#     ticks = np.arange(0, 162, 10)
#     ax.set_xticks(ticks)
#     ax.set_xticklabels(df.dem_id.to_numpy()[ticks])

#     cbar = fig.colorbar(im, ax=ax, pad=0.01)
#     cbar.set_label("Symmetric DEM discrepancy")

#     finish(fig, "04_python_idl_agreement_heatmap.png")


# def plot_chisq_vs_recovery(df):
#     """Investigate intensity-fit quality versus DEM recovery."""
#     fig, axes = plt.subplots(1, 2, figsize=(13, 5.5), sharey=True)

#     configurations = [
#         ("python_raw_chisq", "python_vs_truth",
#          "Python", "#2878b5"),
#         ("idl_raw_chisq", "idl_vs_truth",
#          "IDL", "#e68628"),
#     ]

#     for ax, (chisq_col, pair, name, color) in zip(
#         axes, configurations
#     ):
#         chisq = df[chisq_col].to_numpy()
#         recovery = distance(df, pair)

#         if np.any(chisq < 0) or not np.isfinite(chisq).all():
#             raise ValueError(f"Invalid {name} chi-square values")

#         positive = chisq[chisq > 0]
#         if len(positive) == 0:
#             raise ValueError(f"No positive {name} chi-square values")

#         # Only for plotting exact zero values on logarithmic axes.
#         floor = positive.min() / 10
#         x = np.maximum(chisq, floor)

#         ax.scatter(x, recovery, s=27, alpha=0.65, color=color)
#         ax.set_xscale("log")
#         ax.set_xlabel(f"{name} raw chi-square")
#         ax.set_title(f"{name}: intensity fit vs synthetic recovery")
#         ax.grid(alpha=0.2, which="both")

#     axes[0].set_ylabel("Gray-band discrepancy from synthetic truth")
#     axes[0].set_ylim(-0.05, 2.05)

#     fig.tight_layout()
#     finish(fig, "05_chisq_vs_truth_recovery.png")


# def load_dem_curves(dem_id):
#     """Load the three actual base-solution DEM curves."""
#     truth = np.loadtxt(
#         DATA / "synthetic_dems_data" / f"DEM_{dem_id}.txt"
#     )

#     idl = readsav(
#         str(DATA / "all_IDL_synthetic_dem_data"
#             / f"idl_synthetic_dem_idx{dem_id}_allfilters_base.sav"),
#         python_dict=True
#     )

#     with np.load(
#         DATA / "output_data" / "xrtpy_base_dem_results"
#         / f"xrtpy_dem_idx{dem_id}_allfilters_base.npz"
#     ) as py:
#         py_logt = np.asarray(py["logT"], dtype=float).ravel()
#         py_dem = np.asarray(py["dem"], dtype=float).ravel()

#     logt = truth[:, 0]
#     synthetic = truth[:, 1]
#     idl_logt = np.asarray(idl["logt_out"], dtype=float).ravel()
#     idl_dem = np.asarray(idl["dem_out"], dtype=float).ravel()

#     for grid in (py_logt, idl_logt):
#         if grid.shape != logt.shape or not np.allclose(
#             grid, logt, rtol=0, atol=1e-5
#         ):
#             raise ValueError(f"DEM {dem_id}: grid mismatch")

#     return logt, synthetic, py_dem, idl_dem


# def plot_actual_dem_curves(df):
#     """Overlay real recovered curves, including difficult cases."""
#     discrepancy = distance(df, "python_vs_idl")
#     ranked = np.argsort(discrepancy)[::-1]

#     selected = [0, 14]
#     for index in ranked:
#         dem_id = int(df.iloc[index]["dem_id"])
#         if dem_id not in selected:
#             selected.append(dem_id)
#         if len(selected) == 4:
#             break

#     fig, axes = plt.subplots(2, 2, figsize=(13, 9))
#     axes = axes.ravel()

#     for ax, dem_id in zip(axes, selected):
#         logt, synthetic, py_dem, idl_dem = load_dem_curves(dem_id)

#         ax.plot(
#             logt, synthetic, "k-", linewidth=2.2,
#             label="Synthetic truth"
#         )
#         ax.plot(
#             logt, py_dem, color="#2878b5",
#             linewidth=1.8, label="Python"
#         )
#         ax.plot(
#             logt, idl_dem, color="#e68628",
#             linewidth=1.8, label="IDL"
#         )

#         # DEM values span many orders of magnitude.
#         ax.set_yscale("log")
#         ax.set_xlabel("log10(T / K)")
#         ax.set_ylabel("DEM (cm^-5 K^-1)")
#         ax.set_title(f"DEM {dem_id}")
#         ax.grid(alpha=0.2)

#     axes[0].legend()
#     fig.suptitle("Recovered DEM curves compared with synthetic truth")
#     fig.tight_layout()

#     finish(fig, "06_actual_dem_curve_examples.png")


# def report_trends(df):
#     """Print interpretable descriptive diagnostics."""
#     py = distance(df, "python_vs_truth")
#     idl = distance(df, "idl_vs_truth")
#     agreement = distance(df, "python_vs_idl")

#     print("\nGRAY-BAND VALIDATION")
#     print("=" * 55)

#     for pair, label, _ in PAIRS:
#         x = distance(df, pair)
#         print(
#             f"{label:25s} "
#             f"median={np.median(x):.4f}, "
#             f"p90={np.percentile(x, 90):.4f}"
#         )

#     print("\nPaired Python vs IDL recovery:")
#     print("Python closer:", int(np.sum(py < idl - 1e-10)))
#     print("IDL closer:", int(np.sum(idl < py - 1e-10)))
#     print("Near-maximum Python distances:",
#           int(np.sum(py >= 2 - 1e-6)))

#     print("\nLargest Python-IDL disagreements:")
#     for idx in np.argsort(agreement)[-10:][::-1]:
#         print(
#             f"DEM {int(df.iloc[idx].dem_id):3d}: "
#             f"Python-IDL={agreement[idx]:.5f}, "
#             f"Python-truth={py[idx]:.5f}, "
#             f"IDL-truth={idl[idx]:.5f}"
#         )

#     print("\nNote: DEM index is an input identifier. Apparent")
#     print("clusters require interpretation of the synthetic")
#     print("DEM families before claiming a physical trend.")


# def main():
#     PLOTS.mkdir(parents=True, exist_ok=True)
#     df = load_metrics()

#     plot_simple_histograms(df)
#     plot_python_vs_idl_recovery(df)
#     plot_case_trends(df)
#     plot_agreement_heatmap(df)
#     plot_chisq_vs_recovery(df)
#     plot_actual_dem_curves(df)
#     report_trends(df)

#     print(f"\nAll diagnostic figures saved to: {PLOTS}")


# if __name__ == "__main__":
#     main()
"""Clear diagnostic plots for three-way synthetic DEM validation.

Run from any working directory:
    python plot_synthetic_triplet_diagnostics.py

Reads existing 162-case validation CSVs; never reruns a DEM inversion.
"""

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data" / "output_data" / "scientific_validation"
OUT = ROOT / "plots" / "scientific_validation" / "diagnostics"

PAIRS = ["python vs truth", "idl vs truth", "python vs idl"]
LABELS = ["Python vs synthetic truth", "IDL vs synthetic truth", "Python vs IDL"]
COLORS = ["#2677B9", "#E18832", "#777777"]
REGIONS = ["Gray band", "Full grid"]


def load_data():
    long = pd.read_csv(DATA / "synthetic_triplet_readable.csv")
    wide = pd.read_csv(DATA / "synthetic_triplet_metrics.csv")
    ids = set(range(162))
    if len(wide) != 162 or set(wide.dem_id) != ids or wide.dem_id.duplicated().any():
        raise ValueError("Expected exactly 162 unique DEM IDs in the wide table")
    expected = {(i, p, r) for i in ids for p in PAIRS for r in REGIONS}
    found = set(zip(long.dem_id, long.comparison, long.region))
    if len(long) != 972 or len(found) != 972 or found != expected:
        raise ValueError("Expected all 972 DEM/pair/region records in the long table")
    for col in ("dem_distance_0_to_2", "em_bias_percent"):
        if not np.isfinite(long[col].to_numpy(dtype=float)).all():
            raise ValueError(f"Non-finite {col}")
    d = long.dem_distance_0_to_2.to_numpy(dtype=float)
    if np.any(d < -1e-10) or np.any(d > 2 + 1e-10):
        raise ValueError("DEM discrepancy outside [0, 2]")
    for col in ("python_raw_chisq", "idl_raw_chisq"):
        a = wide[col].to_numpy(dtype=float)
        if not np.isfinite(a).all() or np.any(a < 0):
            raise ValueError(f"Invalid {col}")
    return long, wide.set_index("dem_id").sort_index()


def series(long, pair, region="Gray band"):
    return (long.loc[(long.comparison == pair) & (long.region == region)]
            .sort_values("dem_id").dem_distance_0_to_2.to_numpy(dtype=float))


def save(fig, filename):
    fig.savefig(OUT / filename, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print("Saved:", OUT / filename)


def histogram(long):
    fig, ax = plt.subplots(figsize=(9, 5.5))
    bins = np.linspace(0, 2.000001, 17)
    for pair, label, color in zip(PAIRS, LABELS, COLORS):
        ax.hist(series(long, pair), bins=bins, alpha=0.44,
                label=label, color=color, edgecolor="white", linewidth=0.45)
    ax.set(xlabel="DEM discrepancy (0 = identical; 2 = maximum)",
           ylabel="Number of synthetic DEMs", title="Gray band: how different are the recovered DEMs?")
    ax.set_xlim(0, 2.02)
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    save(fig, "01_overlapping_histograms_gray_band.png")


def median_bars(long):
    medians = np.array([[np.median(series(long, pair, region)) for pair in PAIRS]
                        for region in REGIONS])
    x = np.arange(3)
    fig, ax = plt.subplots(figsize=(10, 5.5))
    width = 0.35
    for j, (region, shade) in enumerate(zip(REGIONS, [0.85, 0.45])):
        positions = x + (j - 0.5) * width
        bars = ax.bar(positions, medians[j], width, alpha=shade,
                      color=COLORS, edgecolor="black", linewidth=0.45,
                      label=region)
        for b, value in zip(bars, medians[j]):
            ax.text(b.get_x() + b.get_width() / 2, b.get_height() + 0.028,
                    f"{value:.2f}", ha="center", fontsize=9)
    ax.set_xticks(x, LABELS)
    ax.set_ylim(0, max(1.08, medians.max() + 0.15))
    ax.set_ylabel("Median DEM discrepancy (lower is better)")
    ax.set_title("Typical DEM discrepancy: gray band vs full grid")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.2)
    fig.tight_layout()
    save(fig, "02_median_discrepancy_bars.png")


def truth_scatter(long):
    py = series(long, "python vs truth")
    idl = series(long, "idl vs truth")
    better_py = py < idl - 1e-10
    better_idl = idl < py - 1e-10
    fig, ax = plt.subplots(figsize=(7.5, 7))
    ax.scatter(idl[better_py], py[better_py], s=32, alpha=0.72,
               color=COLORS[0], label=f"Python closer: {better_py.sum()}")
    ax.scatter(idl[better_idl], py[better_idl], s=32, alpha=0.72,
               color=COLORS[1], label=f"IDL closer: {better_idl.sum()}")
    ties = ~(better_py | better_idl)
    if ties.any():
        ax.scatter(idl[ties], py[ties], s=32, color=COLORS[2], label=f"Ties: {ties.sum()}")
    ax.plot([0, 2], [0, 2], "k--", linewidth=1)
    ax.set(xlim=(0, 2.04), ylim=(0, 2.04),
           xlabel="IDL vs synthetic truth: DEM discrepancy",
           ylabel="Python vs synthetic truth: DEM discrepancy",
           title="Gray band: which solver is closer to synthetic truth?")
    ax.text(0.04, 0.95, "Above diagonal: IDL closer", transform=ax.transAxes,
            va="top", fontsize=9)
    ax.text(0.96, 0.05, "Below diagonal: Python closer", transform=ax.transAxes,
            ha="right", fontsize=9)
    ax.grid(alpha=0.2)
    ax.legend(frameon=False, loc="center right")
    fig.tight_layout()
    save(fig, "03_truth_recovery_scatter.png")


def paired_differences(long):
    py = series(long, "python vs truth")
    idl = series(long, "idl vs truth")
    diff = py - idl
    ids = np.arange(162)
    fig, ax = plt.subplots(figsize=(12, 5.5))
    ax.bar(ids, diff, width=0.9,
           color=np.where(diff < -1e-10, COLORS[0], COLORS[1]))
    ax.axhline(0, color="black", linewidth=1)
    ax.axhline(np.median(diff), color="black", linestyle="--", alpha=0.6,
               label=f"Median difference: {np.median(diff):+.3f}")
    ax.set(xlim=(-1, 162), xlabel="Synthetic DEM ID",
           ylabel="Python distance minus IDL distance",
           title="Gray band: which solver does better for each synthetic DEM?")
    ax.text(0.01, 0.96, "Above zero: IDL closer", transform=ax.transAxes, va="top")
    ax.text(0.01, 0.04, "Below zero: Python closer", transform=ax.transAxes, va="bottom")
    ax.legend(frameon=False)
    ax.grid(axis="y", alpha=0.18)
    fig.tight_layout()
    save(fig, "04_case_by_case_advantage.png")


def python_idl_vs_chisq(long, wide):
    agreement = series(long, "python vs idl")
    # x is signed log-ratio: > 0 means Python has higher raw chi-square.
    a = wide.python_raw_chisq.to_numpy(dtype=float)
    b = wide.idl_raw_chisq.to_numpy(dtype=float)
    positive = np.r_[a[a > 0], b[b > 0]]
    floor = positive.min() * 0.1
    ratio = np.log10(np.maximum(a, floor) / np.maximum(b, floor))
    fig, ax = plt.subplots(figsize=(9, 6))
    sc = ax.scatter(ratio, agreement, c=np.arange(162), cmap="viridis",
                    s=40, alpha=0.8, edgecolor="none")
    ax.axvline(0, color="black", linestyle="--", linewidth=1)
    ax.set(xlabel="log10(Python raw chi-square / IDL raw chi-square)",
           ylabel="Python vs IDL DEM discrepancy (gray band)",
           ylim=(-0.02, 2.05),
           title="Does intensity-fit disagreement track DEM disagreement?")
    ax.text(0.02, 0.04, "Left: Python lower χ² | Right: IDL lower χ²",
            transform=ax.transAxes, fontsize=9)
    fig.colorbar(sc, ax=ax, label="Synthetic DEM ID")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save(fig, "05_python_idl_distance_vs_chisq_ratio.png")


def python_idl_em_bias(long):
    # Separate signed EM bias, which can be extremely large, from bounded distance.
    subset = (long.loc[(long.comparison == "python vs idl") &
                       (long.region == "Gray band")]
              .sort_values("dem_id"))
    x = subset.dem_distance_0_to_2.to_numpy(dtype=float)
    bias = subset.em_bias_percent.to_numpy(dtype=float)
    # asinh displays large signed percentages without discarding or clipping them.
    transform = np.arcsinh(bias / 100.0)
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.scatter(x, transform, s=36, color=COLORS[0], alpha=0.7)
    ax.axhline(0, color="black", linestyle="--", linewidth=1)
    ticks_pct = np.array([-10000, -1000, -100, -25, 0, 25, 100, 1000, 10000])
    ax.set_yticks(np.arcsinh(ticks_pct / 100.0))
    ax.set_yticklabels([f"{v:+,d}%" if v else "0%" for v in ticks_pct])
    ax.set(xlabel="Python vs IDL DEM discrepancy (gray band)",
           ylabel="Python EM relative to IDL (signed percent; nonlinear scale)",
           title="Python–IDL: DEM shape disagreement versus total EM bias")
    ax.grid(alpha=0.2)
    fig.tight_layout()
    save(fig, "06_python_idl_em_bias_vs_distance.png")


def main():
    long, wide = load_data()
    OUT.mkdir(parents=True, exist_ok=True)
    histogram(long)
    median_bars(long)
    truth_scatter(long)
    paired_differences(long)
    python_idl_vs_chisq(long, wide)
    python_idl_em_bias(long)
    py = series(long, "python vs truth")
    idl = series(long, "idl vs truth")
    print("\nGRAY BAND — 162 CASES")
    print(f"Python closer to synthetic truth: {np.count_nonzero(py < idl - 1e-10)}")
    print(f"IDL closer to synthetic truth: {np.count_nonzero(idl < py - 1e-10)}")
    print(f"Numerical ties: {np.count_nonzero(np.abs(py - idl) <= 1e-10)}")
    print(f"Median Python–IDL discrepancy: {np.median(series(long, 'python vs idl')):.4f}")
    print("\nThe plots show associations, not verified causes of solver behavior.")


if __name__ == "__main__":
    main()

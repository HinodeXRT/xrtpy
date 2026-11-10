
"""Summarize and plot the 162-case synthetic DEM validation.

Reads previously calculated metrics. Does not invoke the DEM solver.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE_DIR = Path(__file__).resolve().parent
RESULTS_DIR = BASE_DIR / "data" / "output_data" / "scientific_validation"
PLOT_DIR = BASE_DIR / "plots" / "scientific_validation"

READABLE_CSV = RESULTS_DIR / "synthetic_triplet_readable.csv"
METRICS_CSV = RESULTS_DIR / "synthetic_triplet_metrics.csv"

COMPARISONS = [
    "python vs truth",
    "idl vs truth",
    "python vs idl",
]
LABELS = [
    "Python vs truth",
    "IDL vs truth",
    "Python vs IDL",
]
COLORS = ["steelblue", "darkorange", "gray"]
REGIONS = ["Gray band", "Full grid"]


def load_results():
    """Load and verify the complete validation result tables."""
    readable = pd.read_csv(READABLE_CSV)
    metrics = pd.read_csv(METRICS_CSV)

    expected_ids = set(range(162))
    if set(metrics["dem_id"]) != expected_ids or len(metrics) != 162:
        raise ValueError("Expected exactly 162 unique DEM cases.")

    if len(readable) != 162 * 3 * 2:
        raise ValueError("Expected exactly 972 comparison rows.")

    keys = ["dem_id", "comparison", "region"]
    if readable.duplicated(keys).any():
        raise ValueError("Duplicate comparison records found.")

    for dem_id in expected_ids:
        subset = readable.loc[readable["dem_id"] == dem_id]
        actual = set(zip(subset["comparison"], subset["region"]))
        expected = {
            (comparison, region)
            for comparison in COMPARISONS
            for region in REGIONS
        }
        if actual != expected:
            raise ValueError(f"DEM {dem_id}: incomplete comparisons.")

    numeric_columns = [
        "dem_distance_0_to_2",
        "em_bias_percent",
    ]
    if not np.isfinite(readable[numeric_columns].to_numpy()).all():
        raise ValueError("Nonfinite validation metrics.")

    distances = readable["dem_distance_0_to_2"].to_numpy()
    if np.any(distances < 0) or np.any(distances > 2 + 1e-12):
        raise ValueError("DEM distances outside their valid range.")

    if not np.isfinite(
        metrics[["python_raw_chisq", "idl_raw_chisq"]].to_numpy()
    ).all():
        raise ValueError("Nonfinite chi-square values.")

    return readable, metrics


def make_summary(readable):
    """Calculate case-level descriptive statistics."""
    rows = []

    for region in REGIONS:
        for comparison in COMPARISONS:
            subset = readable.loc[
                (readable["region"] == region)
                & (readable["comparison"] == comparison)
            ]
            distances = subset["dem_distance_0_to_2"].to_numpy()
            biases = subset["em_bias_percent"].to_numpy()

            rows.append({
                "region": region,
                "comparison": comparison,
                "n_cases": len(subset),
                "distance_median": np.median(distances),
                "distance_p25": np.percentile(distances, 25),
                "distance_p75": np.percentile(distances, 75),
                "distance_p90": np.percentile(distances, 90),
                "em_bias_median_percent": np.median(biases),
                "abs_em_bias_median_percent": np.median(np.abs(biases)),
            })

    summary = pd.DataFrame(rows)
    path = RESULTS_DIR / "synthetic_triplet_summary.csv"
    summary.to_csv(path, index=False, float_format="%.8g")
    print(f"Saved: {path}")

    print("\nDEM discrepancy summary (lower is better):")
    print(summary[
        ["region", "comparison", "distance_median",
         "distance_p25", "distance_p75", "distance_p90"]
    ].to_string(index=False, float_format=lambda x: f"{x:.4f}"))

    return summary


def plot_gray_band_distances(readable):
    """Plot all three gray-band discrepancy distributions."""
    data = [
        readable.loc[
            (readable["region"] == "Gray band")
            & (readable["comparison"] == comparison),
            "dem_distance_0_to_2",
        ].to_numpy()
        for comparison in COMPARISONS
    ]

    fig, ax = plt.subplots(figsize=(9, 6))
    parts = ax.boxplot(
        data, tick_labels=LABELS, patch_artist=True,
        showfliers=True, whis=(5, 95),
    )

    for patch, color in zip(parts["boxes"], COLORS):
        patch.set_facecolor(color)
        patch.set_alpha(0.5)

    ax.set_ylabel("Symmetric DEM discrepancy (0 = identical, 2 = maximum)")
    ax.set_title("Synthetic validation: gray contribution band")
    ax.set_ylim(bottom=0, top=2.05)
    ax.grid(axis="y", alpha=0.25)

    fig.tight_layout()
    fig.savefig(
        PLOT_DIR / "gray_band_dem_discrepancy.png", dpi=200
    )
    plt.close(fig)


def plot_band_vs_full(readable):
    """Compare each case's gray-band and full-grid discrepancies."""
    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharex=True, sharey=True)

    for ax, comparison, label, color in zip(
        axes, COMPARISONS, LABELS, COLORS
    ):
        subset = readable.loc[
            readable["comparison"] == comparison
        ]
        pivot = subset.pivot(
            index="dem_id",
            columns="region",
            values="dem_distance_0_to_2",
        )

        ax.scatter(
            pivot["Full grid"], pivot["Gray band"],
            s=19, alpha=0.65, color=color,
        )
        ax.plot([0, 2], [0, 2], "k--", linewidth=1)
        ax.set_title(label)
        ax.set_xlim(0, 2.05)
        ax.set_ylim(0, 2.05)
        ax.set_xlabel("Full-grid discrepancy")
        ax.grid(alpha=0.2)

    axes[0].set_ylabel("Gray-band discrepancy")
    fig.suptitle("Contribution band versus full temperature grid")
    fig.tight_layout()
    fig.savefig(
        PLOT_DIR / "gray_band_vs_full_grid.png", dpi=200
    )
    plt.close(fig)


def plot_chisq(metrics):
    """Show paired raw chi-square values on logarithmic axes."""
    py = metrics["python_raw_chisq"].to_numpy()
    idl = metrics["idl_raw_chisq"].to_numpy()

    if np.any(py < 0) or np.any(idl < 0):
        raise ValueError("Raw chi-square cannot be negative.")

    positive = np.concatenate([py[py > 0], idl[idl > 0]])
    if positive.size == 0:
        raise ValueError("No positive chi-square values to plot.")

    # Plotting floor only, not used for scientific calculations.
    floor = float(positive.min()) / 10.0
    x = np.maximum(idl, floor)
    y = np.maximum(py, floor)

    lower = min(x.min(), y.min()) / 2.0
    upper = max(x.max(), y.max()) * 2.0

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.scatter(x, y, s=24, alpha=0.7, color="steelblue")
    ax.plot(
        [lower, upper], [lower, upper],
        "k--", linewidth=1, label="Equal chi-square",
    )

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(lower, upper)
    ax.set_ylim(lower, upper)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("IDL raw chi-square")
    ax.set_ylabel("Python raw chi-square")
    ax.set_title("Base-solution intensity-fit comparison")
    ax.legend()
    ax.grid(alpha=0.2, which="both")

    fig.tight_layout()
    fig.savefig(
        PLOT_DIR / "python_vs_idl_chisq.png", dpi=200
    )
    plt.close(fig)


def main():
    PLOT_DIR.mkdir(parents=True, exist_ok=True)

    readable, metrics = load_results()
    make_summary(readable)

    plot_gray_band_distances(readable)
    plot_band_vs_full(readable)
    plot_chisq(metrics)

    print(f"\nSaved three figures to: {PLOT_DIR}")
    print("Scientific summary completed.")


if __name__ == "__main__":
    main()

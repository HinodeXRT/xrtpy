"""Reproduce the Part 4 and Part 5 DEM validation summary figures.

Run in synthetic_dem_validation/ after generating synthetic_triplet_metrics.csv.
No inversion runs are performed.
"""
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent
# Supports running a downloaded script placed in synthetic_dem_validation/.
INPUT = ROOT / 'data/output_data/scientific_validation/synthetic_triplet_metrics.csv'
OUTPUT = ROOT / 'plots/scientific_validation'
PAIRS = ['python_vs_truth', 'idl_vs_truth', 'python_vs_idl']
LABELS = ['Python vs truth', 'IDL vs truth', 'Python vs IDL']
REGIONS = [('gray_band', 'Gray band'), ('full_grid', 'Full grid')]


def main():
    df = pd.read_csv(INPUT).sort_values('dem_id')
    if len(df) != 162 or df.dem_id.nunique() != 162 or set(df.dem_id) != set(range(162)):
        raise ValueError('Expected exactly one row for each DEM ID 0..161')
    OUTPUT.mkdir(parents=True, exist_ok=True)
    medians = np.array([[np.median(df[f'{p}_{r}_distance']) for p in PAIRS]
                        for r, _ in REGIONS])

    # Figure 1: Part 4, primary gray-band medians.
    fig, ax = plt.subplots(figsize=(8, 5))
    bars = ax.bar(LABELS, medians[0])
    ax.bar_label(bars, fmt='%.4f', padding=4)
    ax.set_ylim(0, max(medians[0]) * 1.23)
    ax.set_ylabel('Median symmetric DEM discrepancy (lower = closer)')
    ax.set_title('162 synthetic DEMs: gray-band median discrepancy')
    ax.grid(axis='y', alpha=0.25)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(OUTPUT / 'part4_gray_band_medians.png', dpi=220)
    plt.close(fig)

    # Figure 2: Part 4, gray vs full grid medians.
    fig, ax = plt.subplots(figsize=(9, 5))
    x = np.arange(len(PAIRS))
    width = 0.36
    b1 = ax.bar(x - width/2, medians[0], width, label='Gray band')
    b2 = ax.bar(x + width/2, medians[1], width, label='Full grid')
    ax.bar_label(b1, fmt='%.3f', padding=3)
    ax.bar_label(b2, fmt='%.3f', padding=3)
    ax.set_xticks(x, LABELS)
    ax.set_ylim(0, max(medians.ravel()) * 1.25)
    ax.set_ylabel('Median symmetric DEM discrepancy')
    ax.set_title('Gray band and full grid: 162-case medians')
    ax.legend()
    ax.grid(axis='y', alpha=0.25)
    ax.set_axisbelow(True)
    fig.tight_layout()
    fig.savefig(OUTPUT / 'part4_gray_vs_full_medians.png', dpi=220)
    plt.close(fig)

    # Figure 3: Part 5, counts from direct paired comparison.
    counts = []
    for region, _ in REGIONS:
        py = df[f'python_vs_truth_{region}_distance'].to_numpy()
        idl = df[f'idl_vs_truth_{region}_distance'].to_numpy()
        diff = py - idl
        py_wins = int(np.count_nonzero(diff < -1e-10))
        idl_wins = int(np.count_nonzero(diff > 1e-10))
        ties = int(np.count_nonzero(np.abs(diff) <= 1e-10))
        counts.append((py_wins, idl_wins, ties))
        print(f'{region}: Python closer={py_wins}, IDL closer={idl_wins}, ties={ties}, '
              f'median paired delta={np.median(diff):.4f}')
    counts = np.array(counts)
    fig, ax = plt.subplots(figsize=(8, 4.5))
    y = np.arange(2)
    py_bar = ax.barh(y, counts[:, 0], label='Python closer to truth')
    idl_bar = ax.barh(y, counts[:, 1], left=counts[:, 0], label='IDL closer to truth')
    if counts[:, 2].any():
        ax.barh(y, counts[:, 2], left=counts[:, 0]+counts[:, 1], label='Tie')
    for idx, (p, i, t) in enumerate(counts):
        if p: ax.text(p/2, idx, str(p), ha='center', va='center')
        if i: ax.text(p+i/2, idx, str(i), ha='center', va='center')
    ax.set_yticks(y, [label for _, label in REGIONS])
    ax.invert_yaxis()
    ax.set_xlim(0, 162)
    ax.set_xlabel('Number of synthetic DEM cases (of 162)')
    ax.set_title('Which reconstruction is closer to truth?')
    ax.legend(loc='upper center', bbox_to_anchor=(0.5, -0.17), ncol=2)
    fig.tight_layout()
    fig.savefig(OUTPUT / 'part5_paired_wins.png', dpi=220, bbox_inches='tight')
    plt.close(fig)
    print(f'Saved 3 charts to {OUTPUT}')


if __name__ == '__main__':
    main()

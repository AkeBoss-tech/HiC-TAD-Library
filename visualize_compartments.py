"""Visualise A/B compartments from Mouse Micro-C data.

Produces:
    media/compartments_sox11_chr12.png
    media/compartments_mir9_chr13.png
    media/compartments_saddle_chr12.png

Uses the real mcool file at 25 kb resolution (compartments are poorly
defined at 5 kb due to sparse coverage per bin).
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cooler

from src.compartments import compute_compartments, compute_oe_matrix, compartment_saddle


FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
RESOLUTION = 25_000   # 25 kb — standard for compartment analysis

REGIONS = {
    'sox11_chr12': 'chr12:25,000,000-30,000,000',
    'mir9_chr13':  'chr13:83,000,000-85,000,000',
}

os.makedirs('media', exist_ok=True)


def plot_compartments(clr, region, label, out_path):
    print(f"  Computing compartments for {label} …")
    try:
        comp_df = compute_compartments(clr, region)
        oe, _ = compute_oe_matrix(clr, region)
    except Exception as exc:
        print(f"  ERROR computing compartments for {label}: {exc}")
        return

    e1 = comp_df['E1'].values
    chrom = comp_df['chrom'].iloc[0]
    starts = comp_df['start'].values / 1e6   # Mb

    fig, axes = plt.subplots(
        3, 1, figsize=(12, 10),
        gridspec_kw={'height_ratios': [0.08, 0.5, 0.5]},
    )

    # ---- Top: compartment bar ----
    ax_bar = axes[0]
    compartment_colors = ['#d62728' if c == 'A' else '#1f77b4' if c == 'B' else '#aaaaaa'
                          for c in comp_df['compartment']]
    ax_bar.bar(
        np.arange(len(e1)), 1,
        color=compartment_colors, width=1, edgecolor='none',
    )
    ax_bar.set_xlim(0, len(e1))
    ax_bar.set_yticks([])
    ax_bar.set_xticks([])
    ax_bar.set_ylabel('Comp.', fontsize=8)

    # ---- Middle: O/E heatmap ----
    ax_oe = axes[1]
    vmax = np.nanpercentile(oe[oe > 0], 99) if (oe > 0).any() else 2.0
    im = ax_oe.imshow(
        oe, aspect='auto', origin='upper',
        cmap='RdBu_r', vmin=0, vmax=vmax,
        extent=[starts[0], starts[-1], starts[-1], starts[0]],
    )
    plt.colorbar(im, ax=ax_oe, shrink=0.8, label='O/E contact frequency')
    ax_oe.set_ylabel('Genomic position (Mb)')
    ax_oe.set_title(f'Observed/Expected contact matrix — {label.replace("_", " ").title()}')

    # ---- Bottom: E1 eigenvector track ----
    ax_e1 = axes[2]
    pos = (comp_df['start'].values + comp_df['end'].values) / 2 / 1e6
    pos_a = np.where(comp_df['compartment'].values == 'A', e1, np.nan)
    pos_b = np.where(comp_df['compartment'].values == 'B', e1, np.nan)

    ax_e1.fill_between(pos, pos_a, 0, color='#d62728', alpha=0.7, label='A compartment')
    ax_e1.fill_between(pos, pos_b, 0, color='#1f77b4', alpha=0.7, label='B compartment')
    ax_e1.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax_e1.set_xlabel('Genomic position (Mb)')
    ax_e1.set_ylabel('E1 eigenvector')
    ax_e1.legend(loc='upper right', fontsize=9)
    ax_e1.set_xlim(pos[0], pos[-1])

    # Summary stats
    a_frac = (comp_df['compartment'] == 'A').mean()
    b_frac = (comp_df['compartment'] == 'B').mean()
    ax_e1.set_title(
        f'First eigenvector — A: {a_frac:.0%}  B: {b_frac:.0%}  ({RESOLUTION//1000} kb bins)',
        fontsize=9,
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_saddle(clr, region, label, out_path):
    print(f"  Computing saddle plot for {label} …")
    try:
        saddle, quantile_edges = compartment_saddle(clr, region, n_quantiles=20)
    except Exception as exc:
        print(f"  ERROR computing saddle for {label}: {exc}")
        return

    fig, ax = plt.subplots(figsize=(6, 5))
    vmax = np.nanpercentile(np.abs(saddle), 95)
    im = ax.imshow(
        np.log2(np.clip(saddle, 0.1, None)),
        aspect='auto', cmap='RdBu_r',
        vmin=-np.log2(vmax + 1), vmax=np.log2(vmax + 1),
        origin='upper',
    )
    plt.colorbar(im, ax=ax, label='log2(mean O/E)')
    n = saddle.shape[0]
    ax.set_xticks([0, n // 2, n - 1])
    ax.set_xticklabels(['B\n(low E1)', 'mid', 'A\n(high E1)'])
    ax.set_yticks([0, n // 2, n - 1])
    ax.set_yticklabels(['B', 'mid', 'A'])
    ax.set_title(f'Compartment saddle plot — {label.replace("_", " ").title()}')
    ax.set_xlabel('E1 quantile (sorted bins)')
    ax.set_ylabel('E1 quantile (sorted bins)')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    if not os.path.exists(FILE_PATH):
        print(f"ERROR: {FILE_PATH} not found. Please place the mcool file at that path.")
        sys.exit(1)

    print(f"Loading {FILE_PATH} at {RESOLUTION // 1000} kb …")
    clr = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION}")

    for label, region in REGIONS.items():
        plot_compartments(clr, region, label, f"media/compartments_{label}.png")

    # Saddle plot for first region
    label0 = list(REGIONS.keys())[0]
    plot_saddle(clr, list(REGIONS.values())[0], label0, f"media/compartments_saddle_{label0}.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  O/E matrix  : Red = contacts more frequent than expected; blue = depleted.")
    print("  E1 track    : Positive (red) = A compartment — active, gene-dense, accessible.")
    print("                Negative (blue) = B compartment — inactive, gene-poor, heterochromatic.")
    print("  Saddle plot : Strong A–A (top-right) and B–B (bottom-left) enrichment indicates")
    print("                well-separated compartments.  B–A mixing (off-diagonal) reflects")
    print("                compartment boundaries or transitions.")


if __name__ == '__main__':
    main()

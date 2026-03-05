"""Chromatin loop calling and visualisation from Mouse Micro-C data.

Produces:
    media/loops_sox11_chr12.png    — contact matrix + called loops
    media/loops_mir9_chr13.png
    media/loops_pileup_sox11.png   — aggregate loop pileup

Uses 5 kb resolution. Loops are called with a donut background model
and BH-corrected z-score filter.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import cooler

from src.loops import call_loops, loop_pileup


FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
RESOLUTION = 5_000

REGIONS = {
    'sox11_chr12': 'chr12:26,000,000-28,000,000',
    'mir9_chr13':  'chr13:83,500,000-84,500,000',
}

os.makedirs('media', exist_ok=True)


def plot_loops(clr, region, label, out_path):
    from src.tad_boundaries import parse_coordinates
    chrom, start, end = parse_coordinates(region)

    print(f"  Calling loops for {label} …")
    try:
        loops_df = call_loops(
            clr, region,
            min_dist_bp=50_000,
            max_dist_bp=1_500_000,
            enrichment_threshold=1.3,
            zscore_threshold=1.5,
            fdr_threshold=0.2,
        )
    except Exception as exc:
        print(f"  ERROR calling loops for {label}: {exc}")
        return

    print(f"  Found {len(loops_df)} loops in {label}")

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix_log = np.log1p(np.nan_to_num(matrix, nan=0.0))

    n = matrix.shape[0]
    resolution = clr.binsize
    extent_mb = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]

    fig, ax = plt.subplots(figsize=(9, 8))
    im = ax.imshow(
        matrix_log, aspect='equal', origin='upper',
        cmap='YlOrRd', extent=extent_mb,
    )
    plt.colorbar(im, ax=ax, shrink=0.8, label='log(1 + balanced contact)')

    # Overlay loop calls
    for _, row in loops_df.iterrows():
        x = (row['start2'] + row['end2']) / 2 / 1e6
        y = (row['start1'] + row['end1']) / 2 / 1e6
        size = max(20, row['enrichment'] * 30)
        ax.scatter(x, y, s=size, c='#0000ff', alpha=0.7, marker='o',
                   linewidths=1.5, edgecolors='white', zorder=5)
        ax.scatter(y, x, s=size, c='#0000ff', alpha=0.7, marker='o',
                   linewidths=1.5, edgecolors='white', zorder=5)

    ax.set_xlabel('Genomic position (Mb)')
    ax.set_ylabel('Genomic position (Mb)')
    ax.set_title(
        f'Chromatin loops — {label.replace("_", " ").title()}\n'
        f'{len(loops_df)} loops called  ({RESOLUTION // 1000} kb resolution)',
    )

    patch = mpatches.Patch(color='#0000ff', label=f'Loop anchors (n={len(loops_df)})')
    ax.legend(handles=[patch], loc='upper left', fontsize=9)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")
    return loops_df


def plot_pileup(clr, loops_df, label, out_path, flank_bins=15):
    if loops_df is None or loops_df.empty:
        print(f"  No loops to pileup for {label}, skipping.")
        return

    print(f"  Computing loop pileup for {label} …")
    try:
        pileup = loop_pileup(clr, loops_df, flank_bins=flank_bins)
    except Exception as exc:
        print(f"  ERROR computing pileup for {label}: {exc}")
        return

    size = 2 * flank_bins + 1
    resolution = clr.binsize
    ticks = np.linspace(-flank_bins, flank_bins, 5) * resolution / 1000
    tick_pos = np.linspace(0, size - 1, 5)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(
        pileup, aspect='equal', origin='upper', cmap='YlOrRd',
        vmin=np.nanpercentile(pileup, 5),
        vmax=np.nanpercentile(pileup, 99),
    )
    plt.colorbar(im, ax=ax, label='Mean balanced contact')
    ax.set_xticks(tick_pos)
    ax.set_xticklabels([f'{t:+.0f}' for t in ticks], fontsize=8)
    ax.set_yticks(tick_pos)
    ax.set_yticklabels([f'{t:+.0f}' for t in ticks], fontsize=8)
    ax.set_xlabel('Distance from anchor 2 (kb)')
    ax.set_ylabel('Distance from anchor 1 (kb)')
    ax.set_title(
        f'Loop pileup — {label.replace("_", " ").title()}\n'
        f'Averaged over {len(loops_df)} loops ± {flank_bins * resolution // 1000} kb',
    )
    ax.axhline(flank_bins, color='white', linewidth=0.5, linestyle='--', alpha=0.5)
    ax.axvline(flank_bins, color='white', linewidth=0.5, linestyle='--', alpha=0.5)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def main():
    if not os.path.exists(FILE_PATH):
        print(f"ERROR: {FILE_PATH} not found.")
        sys.exit(1)

    print(f"Loading {FILE_PATH} at {RESOLUTION // 1000} kb …")
    clr = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION}")

    all_loops = {}
    for label, region in REGIONS.items():
        loops = plot_loops(clr, region, label, f"media/loops_{label}.png")
        all_loops[label] = loops

    # Pileup for first region
    label0 = list(REGIONS.keys())[0]
    plot_pileup(clr, all_loops.get(label0), label0, f"media/loops_pileup_{label0}.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  Blue dots  : Called loop anchors (size scales with enrichment score).")
    print("               Loops appear as dots OFF the main diagonal — each dot marks")
    print("               two genomic loci in 3D contact (loop base).")
    print("  Pileup     : Aggregate contact matrix centred on each pair of loop anchors.")
    print("               A strong signal at the centre confirms bona-fide point-enrichment")
    print("               (corner peak) characteristic of CTCF/cohesin-mediated loops.")


if __name__ == '__main__':
    main()

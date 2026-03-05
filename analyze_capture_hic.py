"""Capture Hi-C interaction analysis (CHiCAGO-inspired).

Uses the real mcool file for the contact matrix.
Capture baits are loaded from data/raw/capture_baits.bed if present.
Otherwise, synthetic baits are placed at strong TAD boundaries
(simulating a CTCF-targeted capture experiment).

Produces:
    media/capture_interactions_sox11.png  — bait-prey interaction arcs
    media/capture_interactions_mir9.png
    media/capture_summary.png             — top bait interaction counts
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import pandas as pd
import cooler

from src.capture_hic import load_baits, call_bait_interactions, summarise_interactions


FILE_PATH   = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
BAITS_PATH  = "data/raw/capture_baits.bed"
RESOLUTION  = 5_000

REGIONS = {
    'sox11_chr12': 'chr12:26,000,000-28,000,000',
    'mir9_chr13':  'chr13:83,500,000-84,500,000',
}

os.makedirs('media', exist_ok=True)


# ---------------------------------------------------------------------------
# Synthetic baits
# ---------------------------------------------------------------------------

def generate_synthetic_baits(clr, region, n_baits=8, seed=42):
    """Place synthetic capture baits at strong TAD boundaries."""
    import cooltools
    from src.tad_boundaries import parse_coordinates, make_view_df, score_boundary_prominence

    chrom, start, end = parse_coordinates(region)
    view_df = make_view_df(chrom, start, end)
    window_bp = max(5 * clr.binsize, 25_000)

    try:
        ins = cooltools.insulation(clr, [window_bp], view_df=view_df)
        ins = ins[(ins['chrom'] == chrom) & (ins['start'] >= start) & (ins['end'] <= end)]
        scored = score_boundary_prominence(ins, window_bp)
        boundaries = scored[scored['boundary_class'] == 'strong'].head(n_baits)
        if len(boundaries) == 0:
            boundaries = scored[scored['boundary_class'].isin(['strong', 'weak'])].head(n_baits)

        if len(boundaries) > 0:
            baits = boundaries[['chrom', 'start', 'end']].copy()
            mid = (baits['start'] + baits['end']) // 2
            baits['start'] = (mid - 500).clip(lower=start)
            baits['end'] = (mid + 500).clip(upper=end)
            print(f"  Generated {len(baits)} synthetic baits at TAD boundaries ({region})")
            return baits.reset_index(drop=True)
    except Exception:
        pass

    # Fallback: evenly spaced
    rng = np.random.default_rng(seed)
    positions = np.linspace(start + 100_000, end - 100_000, n_baits).astype(int)
    df = pd.DataFrame({
        'chrom': chrom,
        'start': positions - 500,
        'end': positions + 500,
    })
    print(f"  Generated {len(df)} evenly-spaced synthetic baits ({region})")
    return df


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_arc_interactions(clr, baits_df, region, label, out_path):
    from src.tad_boundaries import parse_coordinates
    chrom, start, end = parse_coordinates(region)

    print(f"  Calling bait interactions for {label} …")
    try:
        idf = call_bait_interactions(
            clr, baits_df, region,
            min_dist_bp=clr.binsize * 2,
            max_dist_bp=1_000_000,
            min_contacts=0.005,
            fdr_threshold=0.1,
        )
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return None

    sig = idf[idf['significant']]
    print(f"  Significant interactions: {len(sig)} / {len(idf)} tested  ({label})")

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix_log = np.log1p(np.nan_to_num(matrix, nan=0.0))
    extent = [start / 1e6, end / 1e6, end / 1e6, start / 1e6]

    fig, axes = plt.subplots(2, 1, figsize=(13, 11),
                              gridspec_kw={'height_ratios': [2, 1]})

    # ---- Top: contact matrix + arcs ----
    ax_mat = axes[0]
    im = ax_mat.imshow(matrix_log, aspect='equal', origin='upper',
                       cmap='YlOrRd', extent=extent)
    plt.colorbar(im, ax=ax_mat, shrink=0.7, label='log(1+balanced)')

    # Bait regions
    region_baits = baits_df[
        (baits_df['chrom'] == chrom) &
        (baits_df['start'] >= start) &
        (baits_df['end'] <= end)
    ]
    for _, bait in region_baits.iterrows():
        mid = (bait['start'] + bait['end']) / 2 / 1e6
        ax_mat.axvline(mid, color='blue', linewidth=1.2, alpha=0.6)
        ax_mat.axhline(mid, color='blue', linewidth=1.2, alpha=0.6)

    for _, row in sig.iterrows():
        bx = (row['bait_start'] + row['bait_end']) / 2 / 1e6
        px = (row['prey_start'] + row['prey_end']) / 2 / 1e6
        alpha = min(0.9, max(0.2, row['enrichment'] / 10))
        ax_mat.scatter(px, bx, s=15, c='#d62728', alpha=alpha, zorder=5)
        ax_mat.scatter(bx, px, s=15, c='#d62728', alpha=alpha, zorder=5)

    ax_mat.set_title(
        f'Capture Hi-C interactions — {label.replace("_", " ").title()}\n'
        f'{len(region_baits)} baits  |  {len(sig)} significant interactions',
    )
    ax_mat.set_xlabel('Genomic position (Mb)')
    ax_mat.set_ylabel('Genomic position (Mb)')

    handles = [
        mpatches.Patch(color='blue', alpha=0.6, label='Bait regions'),
        mpatches.Patch(color='#d62728', label=f'Sig. interactions (n={len(sig)})'),
    ]
    ax_mat.legend(handles=handles, loc='upper left', fontsize=9)

    # ---- Bottom: arc plot (viewpoint from first bait) ----
    ax_arc = axes[1]
    first_bait = region_baits.iloc[0] if not region_baits.empty else None

    if first_bait is not None:
        bait_mid = (first_bait['start'] + first_bait['end']) / 2
        bait_sig = sig[
            (sig['bait_start'] == first_bait['start']) &
            (sig['bait_end'] == first_bait['end'])
        ]

        pos_mb = np.arange(start, end, clr.binsize) / 1e6
        matrix_row = matrix[int((bait_mid - start) // clr.binsize), :]
        if len(matrix_row) > len(pos_mb):
            matrix_row = matrix_row[:len(pos_mb)]
        elif len(pos_mb) > len(matrix_row):
            pos_mb = pos_mb[:len(matrix_row)]

        ax_arc.fill_between(pos_mb, matrix_row, alpha=0.3, color='steelblue', label='Contact profile')
        ax_arc.axvline(bait_mid / 1e6, color='blue', linewidth=2, label=f'Bait ({bait_mid/1e6:.2f} Mb)')

        for _, row in bait_sig.iterrows():
            prey_mid = (row['prey_start'] + row['prey_end']) / 2 / 1e6
            ax_arc.axvline(prey_mid, color='#d62728', linewidth=1, alpha=0.7)

        ax_arc.set_xlabel('Genomic position (Mb)')
        ax_arc.set_ylabel('Balanced contact')
        ax_arc.set_title(f'Virtual 4C from bait at {bait_mid/1e6:.2f} Mb  |  {len(bait_sig)} interactions')
        ax_arc.legend(fontsize=8)
    else:
        ax_arc.text(0.5, 0.5, 'No baits in region', ha='center', va='center',
                    transform=ax_arc.transAxes)
        ax_arc.axis('off')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")
    return idf


def plot_summary(all_idf, out_path):
    if not all_idf:
        return

    combined = pd.concat(all_idf, ignore_index=True)
    sig = combined[combined['significant']]

    if sig.empty:
        print("  No significant interactions for summary plot.")
        return

    summary = summarise_interactions(sig, significant_only=False)

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    ax = axes[0]
    top = summary.head(20)
    labels = [f"{r['bait_chrom']}:{r['bait_start']//1e6:.1f}M" for _, r in top.iterrows()]
    ax.barh(range(len(top)), top['n_interactions'].values, color='#2ca02c', edgecolor='none')
    ax.set_yticks(range(len(top)))
    ax.set_yticklabels(labels, fontsize=8)
    ax.set_xlabel('Number of significant interactions')
    ax.set_title('Top baits by interaction count')
    ax.invert_yaxis()

    ax2 = axes[1]
    ax2.scatter(summary['n_interactions'], summary['max_enrichment'],
                color='#9467bd', alpha=0.7, edgecolors='white')
    ax2.set_xlabel('Number of interactions')
    ax2.set_ylabel('Max enrichment score')
    ax2.set_title('Interaction count vs enrichment per bait')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not os.path.exists(FILE_PATH):
        print(f"ERROR: {FILE_PATH} not found.")
        sys.exit(1)

    demo_mode = not os.path.exists(BAITS_PATH)
    if demo_mode:
        print("NOTE: Bait BED not found — generating synthetic baits at TAD boundaries.")
        print(f"  To use real data, place capture bait coordinates at {BAITS_PATH}")
    else:
        print(f"Loading baits from {BAITS_PATH}")

    print(f"Loading {FILE_PATH} at {RESOLUTION // 1000} kb …")
    clr = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION}")

    all_idf = []
    for label, region in REGIONS.items():
        if demo_mode:
            baits_df = generate_synthetic_baits(clr, region)
        else:
            from src.tad_boundaries import parse_coordinates
            chrom, _, _ = parse_coordinates(region)
            baits_df = load_baits(BAITS_PATH, chrom=chrom)

        idf = plot_arc_interactions(clr, baits_df, region, label,
                                    f"media/capture_interactions_{label}.png")
        if idf is not None:
            all_idf.append(idf)

    plot_summary(all_idf, "media/capture_summary.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  Arc/interaction plot : Red dots mark bait-prey pairs whose contact frequency")
    print("    significantly exceeds the distance-matched background.  Each bait (blue line)")
    print("    is the capture probe position; interactions 'radiate' to prey bins.")
    print("  Virtual 4C profile  : Contact profile from one bait across the region.  Peaks")
    print("    in this profile that exceed background = captured interactions; these are")
    print("    the same contacts that form regulatory loops (enhancer–promoter, CTCF–CTCF).")
    print("  Summary plot : Highly-connected baits (many interactions) often sit at gene")
    print("    promoters or major architectural elements like super-enhancers.")


if __name__ == '__main__':
    main()

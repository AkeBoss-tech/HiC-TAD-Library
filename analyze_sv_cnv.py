"""Structural variant and copy number variation analysis.

Hi-C SV detection uses the real mcool file.
BAM-based SV detection uses data/raw/aligned.bam if present; otherwise
synthetic data is generated to demonstrate the pipeline.

Produces:
    media/sv_translocation_heatmap.png   — inter-chrom contact matrix
    media/sv_deletion_scan.png           — deletion candidates in region
    media/cnv_genome_wide.png            — copy number per chromosome
    media/cnv_segments.png               — segmented CN profile for chr12
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import pandas as pd
import cooler

from src.sv_detection import detect_translocations_hic, detect_deletions_hic, detect_inversions_hic
from src.cnv import compute_hic_coverage, normalise_to_ploidy, segment_cnv, call_cnv


FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
BAM_PATH  = "data/raw/aligned.bam"
RESOLUTION = 50_000   # 50 kb for SV/CNV — coarser for genome-wide coverage
RESOLUTION_DETAIL = 10_000  # 10 kb for within-region deletion scan

REGIONS = {
    'sox11_chr12': 'chr12:26,000,000-28,000,000',
    'mir9_chr13':  'chr13:83,500,000-84,500,000',
}
CHROMOSOMES = [f'chr{i}' for i in range(1, 20)]  # mm10 autosomes

os.makedirs('media', exist_ok=True)


# ---------------------------------------------------------------------------
# Translocation heatmap
# ---------------------------------------------------------------------------

def plot_translocation_heatmap(clr, chromosomes, out_path):
    print("  Scanning for translocations (Hi-C inter-chrom contacts) …")
    try:
        sv_df = detect_translocations_hic(
            clr, chromosomes=chromosomes,
            resolution_bp=RESOLUTION,
            zscore_threshold=5.0,
        )
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    print(f"  Found {len(sv_df)} translocation candidates")

    # Build chromosome-level contact counts matrix
    chrom_list = [c for c in chromosomes if c in clr.chromnames]
    n = len(chrom_list)
    count_mat = np.zeros((n, n))

    for _, row in sv_df.iterrows():
        try:
            i = chrom_list.index(row['chrom1'])
            j = chrom_list.index(row['chrom2'])
            count_mat[i, j] += 1
            count_mat[j, i] += 1
        except ValueError:
            pass

    fig, (ax_mat, ax_bar) = plt.subplots(1, 2, figsize=(14, 6),
                                          gridspec_kw={'width_ratios': [2, 1]})

    im = ax_mat.imshow(count_mat, aspect='equal', cmap='Reds', origin='upper')
    plt.colorbar(im, ax=ax_mat, label='# candidate translocation bins')
    ax_mat.set_xticks(range(n))
    ax_mat.set_xticklabels([c.replace('chr', '') for c in chrom_list], fontsize=7)
    ax_mat.set_yticks(range(n))
    ax_mat.set_yticklabels([c.replace('chr', '') for c in chrom_list], fontsize=7)
    ax_mat.set_title('Inter-chromosomal contact enrichment\n(potential translocation hotspots)')

    # Bar chart of top chromosome pairs
    if not sv_df.empty:
        pair_counts = sv_df.groupby(['chrom1', 'chrom2']).size().sort_values(ascending=False).head(10)
        labels = [f"{c1.replace('chr','')}–{c2.replace('chr','')}" for c1, c2 in pair_counts.index]
        ax_bar.barh(range(len(pair_counts)), pair_counts.values, color='#d62728')
        ax_bar.set_yticks(range(len(pair_counts)))
        ax_bar.set_yticklabels(labels, fontsize=9)
        ax_bar.set_xlabel('# candidate bins')
        ax_bar.set_title('Top chromosome pairs\nby translocation signal')
        ax_bar.invert_yaxis()
    else:
        ax_bar.text(0.5, 0.5, 'No translocations\ndetected', ha='center', va='center',
                    transform=ax_bar.transAxes, fontsize=12)
        ax_bar.axis('off')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Deletion scan within a region
# ---------------------------------------------------------------------------

def plot_deletion_scan(clr, region, label, out_path):
    from src.tad_boundaries import parse_coordinates
    chrom, start, end = parse_coordinates(region)

    print(f"  Scanning for deletions in {label} …")

    try:
        del_df = detect_deletions_hic(clr, region)
        inv_df = detect_inversions_hic(clr, region)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    matrix = clr.matrix(balance=False).fetch(f"{chrom}:{start}-{end}").astype(float)
    matrix = np.nan_to_num(matrix, nan=0.0)
    coverage = matrix.sum(axis=1)
    bins_mb = (start + np.arange(len(coverage)) * clr.binsize) / 1e6

    fig, axes = plt.subplots(2, 1, figsize=(12, 7), sharex=True)

    ax_cov = axes[0]
    ax_cov.fill_between(bins_mb, coverage, alpha=0.6, color='steelblue')
    ax_cov.set_ylabel('Row sum (raw contacts)')
    ax_cov.set_title(f'Coverage and deletion scan — {label.replace("_", " ").title()}')

    if not del_df.empty:
        for _, row in del_df.iterrows():
            x0, x1 = row['start'] / 1e6, row['end'] / 1e6
            ax_cov.axvspan(x0, x1, color='red', alpha=0.3, label='Deletion candidate')

    handles = [plt.matplotlib.patches.Patch(color='red', alpha=0.3, label='Deletion candidate')]
    if not del_df.empty:
        ax_cov.legend(handles=handles, fontsize=9)

    ax_inv = axes[1]
    if not inv_df.empty:
        inv_pos_mb = (inv_df['start'].values + inv_df['end'].values) / 2 / 1e6
        ax_inv.vlines(inv_pos_mb, 0, inv_df['flip_score'].values,
                      colors='darkorange', linewidth=2, label='Inversion candidate')
        ax_inv.set_ylabel('DI flip score')
        ax_inv.legend(fontsize=9)
    else:
        ax_inv.text(0.02, 0.5, 'No inversion candidates detected',
                    transform=ax_inv.transAxes, va='center', fontsize=10, color='grey')
        ax_inv.set_ylabel('DI flip score')

    axes[-1].set_xlabel('Genomic position (Mb)')

    stats_txt = (
        f"Deletions: {len(del_df)}  |  Inversion candidates: {len(inv_df)}"
    )
    fig.text(0.5, 0.01, stats_txt, ha='center', fontsize=9, color='#555555')

    plt.tight_layout(rect=[0, 0.03, 1, 1])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}  (deletions={len(del_df)}, inversions={len(inv_df)})")


# ---------------------------------------------------------------------------
# Genome-wide CNV
# ---------------------------------------------------------------------------

def plot_cnv_genome_wide(clr, chromosomes, out_path):
    print("  Computing genome-wide CNV …")
    try:
        seg_df = call_cnv(clr, chromosomes=chromosomes, smooth_bins=10)
    except Exception as exc:
        print(f"  ERROR computing CNV: {exc}")
        return

    print(f"  Found {len(seg_df)} CNV segments across {seg_df['chrom'].nunique()} chromosomes")

    cn_colors = {
        'homozygous_deletion': '#d62728',
        'heterozygous_deletion': '#ff7f0e',
        'diploid': '#2ca02c',
        'gain': '#1f77b4',
        'amplification': '#9467bd',
    }

    fig, axes = plt.subplots(
        len(chromosomes), 1,
        figsize=(16, max(8, len(chromosomes) * 0.5)),
        sharex=False,
    )
    if len(chromosomes) == 1:
        axes = [axes]

    for ax, chrom in zip(axes, chromosomes):
        try:
            cov_df = compute_hic_coverage(clr, chrom)
            cov_df = normalise_to_ploidy(cov_df, smooth_bins=10)
        except Exception:
            ax.set_visible(False)
            continue

        pos_mb = (cov_df['start'].values + cov_df['end'].values) / 2 / 1e6
        cn = cov_df['copy_number'].values
        ax.fill_between(pos_mb, cn, 2, where=cn > 2, color='#1f77b4', alpha=0.5)
        ax.fill_between(pos_mb, cn, 2, where=cn < 2, color='#d62728', alpha=0.5)
        ax.plot(pos_mb, cn, color='black', linewidth=0.5)
        ax.axhline(2, color='grey', linewidth=0.8, linestyle='--', alpha=0.6)
        ax.set_ylim(0, max(4, cn.max() + 0.5))
        ax.set_yticks([0, 2, 4])
        ax.set_yticklabels(['0', '2', '4'], fontsize=6)
        ax.set_ylabel(chrom.replace('chr', ''), fontsize=7, rotation=0, labelpad=20)
        ax.set_xlim(pos_mb[0], pos_mb[-1])

    fig.text(0.5, 0.01, 'Genomic position (Mb)', ha='center', fontsize=11)
    fig.text(0.01, 0.5, 'Copy number', va='center', rotation='vertical', fontsize=11)
    plt.suptitle('Genome-wide copy number variation (Hi-C coverage)\nBlue = gain,  Red = loss,  Dashed = diploid (CN=2)',
                 fontsize=11)
    plt.tight_layout(rect=[0.03, 0.03, 1, 0.97])
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


# ---------------------------------------------------------------------------
# Detailed CN segment plot for one chromosome
# ---------------------------------------------------------------------------

def plot_cnv_segments(clr, chrom, out_path):
    print(f"  Segmenting CNV for {chrom} …")
    try:
        cov_df = compute_hic_coverage(clr, chrom)
        cov_df = normalise_to_ploidy(cov_df, smooth_bins=5)
        seg_df = segment_cnv(cov_df, change_threshold=0.4)
    except Exception as exc:
        print(f"  ERROR: {exc}")
        return

    pos_mb = (cov_df['start'].values + cov_df['end'].values) / 2 / 1e6
    cn = cov_df['copy_number'].values

    cn_colors = {
        'homozygous_deletion': '#d62728',
        'heterozygous_deletion': '#ff7f0e',
        'diploid': '#2ca02c',
        'gain': '#1f77b4',
        'amplification': '#9467bd',
    }

    fig, ax = plt.subplots(figsize=(14, 4))
    ax.scatter(pos_mb, cn, s=2, color='#888888', alpha=0.4, label='Bin CN')

    for _, row in seg_df.iterrows():
        color = cn_colors.get(row['cn_state'], '#555555')
        x0, x1 = row['start'] / 1e6, row['end'] / 1e6
        ax.hlines(row['mean_cn'], x0, x1, colors=color, linewidth=3)

    ax.axhline(2, color='black', linewidth=0.8, linestyle='--', alpha=0.5, label='Diploid (CN=2)')
    ax.set_xlabel('Genomic position (Mb)')
    ax.set_ylabel('Copy number')
    ax.set_title(f'CNV segmentation — {chrom}  ({len(seg_df)} segments)')

    # Legend for CN states
    legend_handles = [
        plt.matplotlib.patches.Patch(color=v, label=k.replace('_', ' ').title())
        for k, v in cn_colors.items()
    ]
    ax.legend(handles=legend_handles, loc='upper right', fontsize=8, ncol=2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}  ({len(seg_df)} segments)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    if not os.path.exists(FILE_PATH):
        print(f"ERROR: {FILE_PATH} not found.")
        sys.exit(1)

    print(f"Loading {FILE_PATH} at {RESOLUTION // 1000} kb …")
    clr = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION}")

    chroms_available = [c for c in CHROMOSOMES if c in clr.chromnames]

    plot_translocation_heatmap(clr, chroms_available, "media/sv_translocation_heatmap.png")

    # Use the higher-res cooler for deletion/inversion scan
    clr_detail = cooler.Cooler(f"{FILE_PATH}::resolutions/{RESOLUTION_DETAIL}")
    for label, region in REGIONS.items():
        plot_deletion_scan(clr_detail, region, label, f"media/sv_deletion_scan_{label}.png")

    plot_cnv_genome_wide(clr, chroms_available, "media/cnv_genome_wide.png")
    plot_cnv_segments(clr, chroms_available[0] if chroms_available else 'chr12', "media/cnv_segments_chr12.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  Translocation heatmap : High-intensity chromosome pairs have more inter-chrom")
    print("    contacts than expected — possible rearrangements or mis-assembly artifacts.")
    print("    True translocations show focal hot-spots; diffuse signal is Hi-C noise.")
    print("  Deletion scan  : Regions with near-zero coverage (gray bars) flanked by elevated")
    print("    bridging contacts (red spans) are candidate deletions.  Coverage drops = missing")
    print("    chromatin material.  In diploid data, only homozygous deletions are visible.")
    print("  CNV profile    : Deviations from CN=2 along a chromosome.  Smooth plateaus above")
    print("    baseline = amplifications; plateaus below = deletions or LOH regions.")


if __name__ == '__main__':
    main()

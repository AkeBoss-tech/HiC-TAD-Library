"""Hi-C based haplotype phasing analysis.

If data/raw/aligned.bam and data/raw/variants.vcf are present, runs the
real Hi-C phasing pipeline.  Otherwise, synthetic het SNPs and read-pair
phasing constraints are generated to demonstrate the analysis.

Produces:
    media/phasing_haplotype_blocks.png  — haplotype block structure
    media/phasing_snp_density.png       — SNP density and phased fraction
    media/phasing_block_lengths.png     — block length distribution
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

from src.tad_boundaries import parse_coordinates

BAM_PATH = "data/raw/aligned.bam"
VCF_PATH = "data/raw/variants.vcf"
REGION   = "chr12:26,000,000-28,000,000"

os.makedirs('media', exist_ok=True)


# ---------------------------------------------------------------------------
# Synthetic demo data
# ---------------------------------------------------------------------------

def generate_synthetic_phasing(region=REGION, n_snps=300, seed=42):
    """Generate synthetic het SNPs and simulated phasing constraints."""
    rng = np.random.default_rng(seed)
    chrom, start, end = parse_coordinates(region)

    # SNPs at random positions
    positions = np.sort(rng.integers(start, end, n_snps))
    bases = list('ACGT')
    refs = rng.choice(bases, n_snps)
    alts = np.array([rng.choice([b for b in bases if b != r]) for r in refs])
    snp_ids = [f"{chrom}:{p}" for p in positions]

    snps_df = pd.DataFrame({
        'chrom': chrom,
        'pos': positions.astype(int),
        'ref': refs,
        'alt': alts,
        'snp_id': snp_ids,
    })

    # Simulate phasing constraints: read pairs linking nearby SNP pairs
    # True haplotype: first half of SNPs on hap0 (then hap1 alternating in blocks)
    block_size = n_snps // 8
    true_hap = np.array([((i // block_size) % 2) for i in range(n_snps)])

    constraints = []
    rng_local = np.random.default_rng(seed + 1)
    for _ in range(n_snps * 5):  # generate many constraints
        i = rng_local.integers(0, n_snps)
        # Pair with a nearby SNP (Hi-C links proximal loci)
        max_j = n_snps - i
        if max_j <= 1:
            continue
        j = i + rng_local.integers(1, min(30, max_j))
        if j >= n_snps:
            continue
        # 90% accuracy (10% error rate simulating noise)
        true_same = (true_hap[i] == true_hap[j])
        observed_same = true_same if rng_local.random() > 0.1 else not true_same
        constraints.append({
            'snp1_id': snp_ids[i],
            'snp1_allele': 'REF' if true_hap[i] == 0 else 'ALT',
            'snp2_id': snp_ids[j],
            'snp2_allele': 'REF' if true_hap[j] == 0 else 'ALT',
            'read_pair_id': f'read_{_}',
            'same_haplotype': observed_same,
        })

    constraints_df = pd.DataFrame(constraints)
    return snps_df, constraints_df, true_hap


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_haplotype_blocks(phased_df, region, out_path):
    chrom, start, end = parse_coordinates(region)
    pos = phased_df['pos'].values / 1e6
    hap = phased_df['haplotype'].values
    block = phased_df['haplotype_block'].values

    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)

    # ---- Haplotype assignment ----
    ax = axes[0]
    hap0 = np.where(hap == 0, 1, np.nan)
    hap1 = np.where(hap == 1, 1, np.nan)
    unphased = np.where(hap == -1, 1, np.nan)

    ax.scatter(pos, hap0, s=8, color='#1f77b4', alpha=0.7, label='Haplotype 0')
    ax.scatter(pos, hap1 * 0, s=8, color='#d62728', alpha=0.7, label='Haplotype 1')
    ax.scatter(pos, np.where(hap == -1, -0.5, np.nan), s=8, color='grey', alpha=0.4, label='Unphased')
    ax.set_yticks([-0.5, 0, 1])
    ax.set_yticklabels(['unphased', 'hap 1', 'hap 0'], fontsize=8)
    ax.set_title('Haplotype assignment per SNP')
    ax.legend(loc='upper right', fontsize=8)

    # ---- Block structure ----
    ax2 = axes[1]
    unique_blocks = [b for b in phased_df['haplotype_block'].unique() if b >= 0]
    cmap = plt.cm.get_cmap('tab20', max(len(unique_blocks), 1))
    for bi, block_id in enumerate(unique_blocks):
        mask = phased_df['haplotype_block'] == block_id
        block_pos = phased_df.loc[mask, 'pos'].values / 1e6
        if len(block_pos) == 0:
            continue
        ax2.axvspan(block_pos[0], block_pos[-1], alpha=0.3, color=cmap(bi % 20))
    ax2.set_ylabel('Block')
    ax2.set_title(f'Haplotype blocks  ({len(unique_blocks)} blocks, '
                  f'{(phased_df["haplotype"] >= 0).mean():.0%} phased)')
    ax2.set_yticks([])

    # ---- SNP density ----
    ax3 = axes[2]
    bin_mb = np.linspace(start / 1e6, end / 1e6, 40)
    density, _ = np.histogram(pos, bins=bin_mb)
    centres = (bin_mb[:-1] + bin_mb[1:]) / 2
    ax3.bar(centres, density, width=(bin_mb[1] - bin_mb[0]) * 0.9,
            color='#2ca02c', edgecolor='none', alpha=0.7)
    ax3.set_xlabel('Genomic position (Mb)')
    ax3.set_ylabel('SNPs per bin')
    ax3.set_title('Het SNP density')

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_block_lengths(phased_df, out_path):
    unique_blocks = [b for b in phased_df['haplotype_block'].unique() if b >= 0]
    if not unique_blocks:
        return

    block_lengths = []
    for block_id in unique_blocks:
        mask = phased_df['haplotype_block'] == block_id
        block_pos = phased_df.loc[mask, 'pos'].values
        if len(block_pos) >= 2:
            block_lengths.append(block_pos[-1] - block_pos[0])

    if not block_lengths:
        return

    block_lengths_kb = np.array(block_lengths) / 1000

    fig, axes = plt.subplots(1, 2, figsize=(11, 5))

    ax = axes[0]
    ax.hist(block_lengths_kb, bins=20, color='#9467bd', edgecolor='white')
    ax.set_xlabel('Block length (kb)')
    ax.set_ylabel('Number of blocks')
    ax.set_title(f'Haplotype block length distribution\n'
                 f'Median: {np.median(block_lengths_kb):.0f} kb  |  '
                 f'N50: {_n50(block_lengths_kb):.0f} kb')
    ax.axvline(np.median(block_lengths_kb), color='red', linestyle='--',
               label=f'Median {np.median(block_lengths_kb):.0f} kb')
    ax.legend(fontsize=9)

    ax2 = axes[1]
    sorted_lengths = np.sort(block_lengths_kb)[::-1]
    cumsum = np.cumsum(sorted_lengths)
    ax2.plot(sorted_lengths, cumsum / cumsum[-1], color='#9467bd', linewidth=2)
    ax2.set_xlabel('Block length (kb)')
    ax2.set_ylabel('Cumulative fraction of phased bases')
    ax2.set_title('Block length CDF')
    ax2.axhline(0.5, color='red', linestyle='--', alpha=0.5, label='N50')
    ax2.legend(fontsize=9)

    plt.suptitle(f'Phasing quality — {len(unique_blocks)} haplotype blocks', fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def plot_snp_density(phased_df, out_path):
    chrom = phased_df['chrom'].iloc[0]
    pos = phased_df['pos'].values / 1e6
    hap = phased_df['haplotype'].values

    phased_frac = (hap >= 0).mean()
    hap0_frac = (hap == 0).mean()
    hap1_frac = (hap == 1).mean()

    fig, ax = plt.subplots(figsize=(8, 5))
    labels = ['Haplotype 0', 'Haplotype 1', 'Unphased']
    sizes  = [hap0_frac, hap1_frac, 1 - phased_frac]
    colors = ['#1f77b4', '#d62728', '#aaaaaa']
    wedges, texts, autotexts = ax.pie(
        sizes, labels=labels, colors=colors,
        autopct='%1.1f%%', startangle=90,
        wedgeprops={'edgecolor': 'white', 'linewidth': 1.5},
    )
    ax.set_title(
        f'Phasing summary — {chrom}:{phased_df["pos"].min()//1e6:.0f}–{phased_df["pos"].max()//1e6:.0f} Mb\n'
        f'{len(phased_df)} het SNPs  |  {phased_frac:.0%} phased',
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {out_path}")


def _n50(lengths):
    """Compute N50 of an array of lengths."""
    s = np.sort(lengths)[::-1]
    cumsum = np.cumsum(s)
    idx = np.searchsorted(cumsum, cumsum[-1] / 2)
    return s[idx] if idx < len(s) else s[-1]


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    demo_mode = not (os.path.exists(BAM_PATH) and os.path.exists(VCF_PATH))

    if demo_mode:
        print("NOTE: BAM/VCF not found — generating synthetic phasing data for demonstration.")
        print(f"  To use real data, place aligned BAM at {BAM_PATH} and VCF at {VCF_PATH}")
        snps_df, constraints_df, _ = generate_synthetic_phasing()
        print(f"  Generated {len(snps_df)} het SNPs, {len(constraints_df)} phasing constraints")

        from src.phasing import phase_snps
        phased_df = phase_snps(constraints_df, snps_df, min_support=2)
    else:
        from src.phasing import run_hic_phasing
        print(f"Running Hi-C phasing for {REGION} …")
        try:
            phased_df, constraints_df = run_hic_phasing(BAM_PATH, VCF_PATH, REGION)
        except Exception as exc:
            print(f"  ERROR: {exc}")
            sys.exit(1)

    n_blocks = phased_df['haplotype_block'].nunique() - (1 if -1 in phased_df['haplotype_block'].values else 0)
    phased_frac = (phased_df['haplotype'] >= 0).mean()

    print(f"  Phasing result: {len(phased_df)} SNPs, {n_blocks} blocks, {phased_frac:.0%} phased")

    plot_haplotype_blocks(phased_df, REGION, "media/phasing_haplotype_blocks.png")
    plot_block_lengths(phased_df, "media/phasing_block_lengths.png")
    plot_snp_density(phased_df, "media/phasing_snp_density.png")

    print("\nDone. Outputs in media/")
    print("\nResult interpretation:")
    print("  Haplotype blocks  : Contiguous SNPs assigned to the same haplotype form a block.")
    print("    Long blocks = many phasing constraints from Hi-C read pairs spanning the region.")
    print("    Blocks break at positions with no informative read pairs (recombination breakpoints")
    print("    or simply low read coverage).")
    print("  Block N50         : Analogous to assembly N50 — half of phased bases lie in blocks")
    print("    longer than this.  Hi-C phasing typically achieves Mb-scale N50 because")
    print("    long-range contacts link distant het SNPs.")
    print("  Unphased SNPs     : SNPs with no linking read pairs remain unphased (-1).")
    print("    These are often in low-complexity regions, centromeres, or segmental duplications.")


if __name__ == '__main__':
    main()

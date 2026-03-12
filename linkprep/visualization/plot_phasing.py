"""
Visualize haplotype phasing results.

Produces:
  phasing_snp_density.png      — SNP density along chromosomes
  phasing_haplotype_blocks.png — block lengths and N50
  phasing_block_lengths.png    — block length distribution
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config
from analysis.phasing import load_vcf_snps, build_phase_blocks, phasing_stats, print_phasing_stats


# ---------------------------------------------------------------------------
# SNP density
# ---------------------------------------------------------------------------

def plot_snp_density(snps_df: pd.DataFrame, out_path: str, bin_size: int = 500_000):
    chroms = snps_df["chrom"].unique()
    fig, axes = plt.subplots(len(chroms), 1, figsize=(12, 3 * len(chroms)),
                              squeeze=False)

    for ax, chrom in zip(axes[:, 0], chroms):
        sub = snps_df[snps_df["chrom"] == chrom]
        chrom_len = sub["pos"].max() + bin_size
        bins_edges = np.arange(0, chrom_len + bin_size, bin_size)
        counts, _ = np.histogram(sub["pos"].values, bins=bins_edges)
        ax.fill_between(bins_edges[:-1] / 1e6, counts, step="post",
                        color="#5C6BC0", alpha=0.8)
        ax.set_ylabel("SNPs/bin")
        ax.set_xlabel(f"Position on {chrom} (Mb)")
        ax.set_title(f"Heterozygous SNP density — {chrom}")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Haplotype block map
# ---------------------------------------------------------------------------

def plot_haplotype_blocks(blocks_df: pd.DataFrame, out_path: str):
    if blocks_df.empty:
        print("  No phase blocks to plot.")
        return

    chroms = blocks_df["chrom"].unique()
    fig, axes = plt.subplots(len(chroms), 1, figsize=(12, 2 * len(chroms) + 1),
                              squeeze=False)

    cmap = plt.cm.tab20
    for ax, chrom in zip(axes[:, 0], chroms):
        sub = blocks_df[blocks_df["chrom"] == chrom].reset_index(drop=True)
        for i, row in sub.iterrows():
            ax.barh(
                0, row["length_bp"] / 1e6,
                left=row["block_start"] / 1e6,
                height=0.6,
                color=cmap(i % 20), edgecolor="white", linewidth=0.3,
            )
        ax.set_yticks([])
        ax.set_xlabel(f"Position on {chrom} (Mb)")
        ax.set_title(f"Haplotype blocks — {chrom}  ({len(sub)} blocks)")

    plt.suptitle("LinkPrep Haplotype Phasing", fontsize=12, y=1.01)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Block length distribution
# ---------------------------------------------------------------------------

def plot_block_lengths(blocks_df: pd.DataFrame, out_path: str):
    if blocks_df.empty:
        return

    lengths_mb = blocks_df["length_bp"].values / 1e6
    stats = phasing_stats(blocks_df)
    n50_mb = stats.get("N50_bp", 0) / 1e6

    fig, axes = plt.subplots(1, 2, figsize=(11, 4))

    # Histogram
    axes[0].hist(lengths_mb, bins=30, color="#7986CB", edgecolor="white")
    axes[0].axvline(n50_mb, color="#E53935", lw=2, ls="--",
                    label=f"N50 = {n50_mb:.1f} Mb")
    axes[0].set_xlabel("Block length (Mb)")
    axes[0].set_ylabel("Count")
    axes[0].set_title("Haplotype block length distribution")
    axes[0].legend()

    # SNPs per block
    axes[1].hist(blocks_df["n_snps"].values, bins=30,
                 color="#4DB6AC", edgecolor="white")
    axes[1].set_xlabel("SNPs per block")
    axes[1].set_ylabel("Count")
    axes[1].set_title("SNPs per haplotype block")

    plt.suptitle("LinkPrep Phasing Statistics", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all():
    os.makedirs(config.FIGURES_DIR, exist_ok=True)

    if not os.path.exists(config.VCF_FILE):
        sys.exit(f"VCF not found: {config.VCF_FILE}\n"
                 "Run: python sample_data/generate_sample_data.py")

    snps_df   = load_vcf_snps(config.VCF_FILE)
    blocks_df = build_phase_blocks(snps_df,
                                   max_gap_bp=500_000,
                                   min_snps=config.PHASING_MIN_SNPS)
    stats     = phasing_stats(blocks_df)
    print_phasing_stats(stats)

    plot_snp_density(
        snps_df,
        os.path.join(config.FIGURES_DIR, "phasing_snp_density.png"),
    )
    plot_haplotype_blocks(
        blocks_df,
        os.path.join(config.FIGURES_DIR, "phasing_haplotype_blocks.png"),
    )
    plot_block_lengths(
        blocks_df,
        os.path.join(config.FIGURES_DIR, "phasing_block_lengths.png"),
    )
    return blocks_df


if __name__ == "__main__":
    run_all()

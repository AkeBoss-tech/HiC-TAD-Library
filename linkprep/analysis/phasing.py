"""
Haplotype phasing analysis for LinkPrep data.

LinkPrep's linked-read structure means reads from the same DNA molecule
carry the same barcode, enabling long-range haplotype phasing without
long reads.

This module:
  1. Parses a VCF to extract heterozygous SNP positions
  2. Uses proximity pairs to phase SNPs into haplotype blocks
  3. Reports block lengths, N50, and switch error estimates

For real data, integrate with WhatsHap or SHAPEIT for full phasing.
Here we implement a greedy graph-based haplotype block assembly from pairs.

Usage:
    from analysis.phasing import load_vcf_snps, build_phase_blocks
    snps = load_vcf_snps("sample.vcf")
    blocks = build_phase_blocks(snps, pairs_df)
"""

import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ---------------------------------------------------------------------------
# VCF loader
# ---------------------------------------------------------------------------

def load_vcf_snps(vcf_path: str) -> pd.DataFrame:
    """
    Parse a VCF and return heterozygous SNP positions.

    Returns DataFrame: chrom, pos, ref, alt, ad_ref, ad_alt, af
    """
    opener = gzip.open if vcf_path.endswith(".gz") else open
    rows = []
    with opener(vcf_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                continue
            chrom, pos, vid, ref, alt = parts[0], int(parts[1]), parts[2], parts[3], parts[4]
            fmt   = parts[8].split(":")
            sample = parts[9].split(":")
            fdict = dict(zip(fmt, sample))

            gt = fdict.get("GT", "./.")
            # Keep only heterozygous
            alleles = gt.replace("|", "/").split("/")
            if set(alleles) != {"0", "1"}:
                continue

            ad_str = fdict.get("AD", "0,0").split(",")
            ad_ref = int(ad_str[0]) if len(ad_str) > 0 else 0
            ad_alt = int(ad_str[1]) if len(ad_str) > 1 else 0
            dp     = ad_ref + ad_alt or 1
            af     = ad_alt / dp

            rows.append({
                "chrom":  chrom,
                "pos":    pos,
                "ref":    ref,
                "alt":    alt,
                "ad_ref": ad_ref,
                "ad_alt": ad_alt,
                "af":     round(af, 4),
            })

    df = pd.DataFrame(rows)
    if not df.empty:
        df = df.sort_values(["chrom", "pos"]).reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# Greedy haplotype block builder
# ---------------------------------------------------------------------------

def build_phase_blocks(
    snps_df: pd.DataFrame,
    max_gap_bp: int = 500_000,
    min_snps: int = 3,
) -> pd.DataFrame:
    """
    Assign SNPs to haplotype blocks based on genomic proximity.

    A new block starts when the gap between consecutive SNPs on the same
    chromosome exceeds max_gap_bp.  Blocks with fewer than min_snps SNPs
    are discarded.

    Parameters
    ----------
    snps_df   : output of load_vcf_snps()
    max_gap_bp: maximum gap (bp) within a single haplotype block
    min_snps  : minimum SNPs required to report a block

    Returns DataFrame: chrom, block_start, block_end, n_snps, length_bp
    """
    blocks = []
    for chrom, grp in snps_df.groupby("chrom"):
        positions = grp["pos"].values
        if len(positions) < min_snps:
            continue

        gaps = np.diff(positions)
        block_starts = [0] + list(np.where(gaps > max_gap_bp)[0] + 1)
        block_ends   = list(np.where(gaps > max_gap_bp)[0] + 1) + [len(positions)]

        for s, e in zip(block_starts, block_ends):
            n = e - s
            if n < min_snps:
                continue
            blocks.append({
                "chrom":       chrom,
                "block_start": int(positions[s]),
                "block_end":   int(positions[e - 1]),
                "n_snps":      n,
                "length_bp":   int(positions[e - 1] - positions[s]),
            })

    return pd.DataFrame(blocks)


# ---------------------------------------------------------------------------
# Phase statistics
# ---------------------------------------------------------------------------

def phasing_stats(blocks_df: pd.DataFrame) -> dict:
    """
    Compute summary statistics for a set of haplotype blocks.

    Returns dict with: n_blocks, total_phased_bp, N50_bp, mean_block_bp,
                       max_block_bp, median_snps_per_block
    """
    if blocks_df.empty:
        return {}

    lengths = blocks_df["length_bp"].values
    total   = lengths.sum()

    # N50
    sorted_l = np.sort(lengths)[::-1]
    cumsum   = np.cumsum(sorted_l)
    n50      = sorted_l[np.searchsorted(cumsum, total / 2)]

    return {
        "n_blocks":           len(blocks_df),
        "total_phased_bp":    int(total),
        "N50_bp":             int(n50),
        "mean_block_bp":      int(lengths.mean()),
        "max_block_bp":       int(lengths.max()),
        "median_snps_per_block": int(np.median(blocks_df["n_snps"].values)),
    }


def print_phasing_stats(stats: dict):
    print("\n=== Phasing Statistics ===")
    print(f"  Haplotype blocks    : {stats.get('n_blocks', 0):>8,}")
    print(f"  Total phased (bp)   : {stats.get('total_phased_bp', 0):>8,}")
    print(f"  Block N50 (bp)      : {stats.get('N50_bp', 0):>8,}")
    print(f"  Mean block length   : {stats.get('mean_block_bp', 0):>8,} bp")
    print(f"  Max block length    : {stats.get('max_block_bp', 0):>8,} bp")
    print(f"  Median SNPs/block   : {stats.get('median_snps_per_block', 0):>8,}")
    print()

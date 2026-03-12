"""
Visualize copy-number variation from LinkPrep coverage data.

Produces:
  cnv_genome_wide.png    — genome-wide log2 ratio plot
  cnv_segments_<chrom>.png — per-chromosome segmentation
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config
from analysis.cnv_analysis import (
    load_coverage, normalize_coverage, call_cnv_segments, summarise_cnv
)


STATE_COLORS = {
    "deletion":      "#D32F2F",
    "loss":          "#FF7043",
    "neutral":       "#90A4AE",
    "gain":          "#43A047",
    "amplification": "#1B5E20",
}


# ---------------------------------------------------------------------------
# Genome-wide plot
# ---------------------------------------------------------------------------

def plot_genome_wide(cov_df: pd.DataFrame, segs_df: pd.DataFrame, out_path: str):
    chroms = cov_df["chrom"].unique()
    n_chroms = len(chroms)
    offsets = {}
    cum = 0
    for chrom in chroms:
        offsets[chrom] = cum
        cum += cov_df[cov_df["chrom"] == chrom]["end"].max()

    fig, ax = plt.subplots(figsize=(14, 4))

    for chrom in chroms:
        sub  = cov_df[cov_df["chrom"] == chrom]
        off  = offsets[chrom]
        xmid = (sub["start"].values + sub["end"].values) / 2 + off
        ax.scatter(xmid, sub["log2_ratio"].values,
                   s=1, alpha=0.3, color="#78909C", rasterized=True)

    for _, row in segs_df.iterrows():
        off  = offsets.get(row["chrom"], 0)
        x1   = row["start"] + off
        x2   = row["end"]   + off
        color = STATE_COLORS.get(row["cn_state"], "#90A4AE")
        ax.plot([x1, x2], [row["mean_log2"], row["mean_log2"]],
                lw=3, color=color, solid_capstyle="butt", alpha=0.9)

    ax.axhline(0,  color="black", lw=0.8, ls="--")
    ax.axhline(1,  color="#43A047", lw=0.5, ls=":", alpha=0.5)
    ax.axhline(-1, color="#D32F2F", lw=0.5, ls=":", alpha=0.5)

    # Chromosome boundary lines
    chrom_mids = []
    for chrom in chroms:
        sub = cov_df[cov_df["chrom"] == chrom]
        off = offsets[chrom]
        boundary = off + sub["end"].max()
        ax.axvline(boundary, color="white", lw=1.5, alpha=0.8)
        chrom_mids.append(off + sub["end"].max() / 2)

    ax.set_xticks(chrom_mids)
    ax.set_xticklabels(list(chroms), fontsize=8, rotation=45)
    ax.set_ylabel("log2(obs/exp)")
    ax.set_title("LinkPrep — Genome-wide CNV Profile", fontsize=12)

    legend_handles = [mpatches.Patch(color=v, label=k)
                      for k, v in STATE_COLORS.items()]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=7,
              ncol=2, framealpha=0.8)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Per-chromosome segmentation
# ---------------------------------------------------------------------------

def plot_chromosome_cnv(cov_df: pd.DataFrame, segs_df: pd.DataFrame,
                         chrom: str, out_path: str):
    sub     = cov_df[cov_df["chrom"] == chrom]
    sub_seg = segs_df[segs_df["chrom"] == chrom]

    xmid = (sub["start"].values + sub["end"].values) / 2

    fig, ax = plt.subplots(figsize=(12, 3))
    ax.scatter(xmid / 1e6, sub["log2_ratio"].values,
               s=4, alpha=0.4, color="#546E7A")

    for _, row in sub_seg.iterrows():
        color = STATE_COLORS.get(row["cn_state"], "#90A4AE")
        ax.plot([row["start"] / 1e6, row["end"] / 1e6],
                [row["mean_log2"], row["mean_log2"]],
                lw=3, color=color, solid_capstyle="butt")

    ax.axhline(0, color="black", lw=0.8, ls="--")
    ax.set_xlabel(f"Position on {chrom} (Mb)")
    ax.set_ylabel("log2(obs/exp)")
    ax.set_title(f"CNV Segments — {chrom}")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all():
    os.makedirs(config.FIGURES_DIR, exist_ok=True)

    if not os.path.exists(config.COVERAGE_BED):
        sys.exit(f"Coverage BED not found: {config.COVERAGE_BED}\n"
                 "Run: python sample_data/generate_sample_data.py")

    cov_df  = load_coverage(config.COVERAGE_BED)
    cov_df  = normalize_coverage(cov_df)
    segs_df = call_cnv_segments(cov_df,
                                smoothing_bins=config.CNV_SMOOTHING_BINS)

    plot_genome_wide(
        cov_df, segs_df,
        os.path.join(config.FIGURES_DIR, "cnv_genome_wide.png"),
    )

    for chrom in cov_df["chrom"].unique():
        plot_chromosome_cnv(
            cov_df, segs_df, chrom,
            os.path.join(config.FIGURES_DIR, f"cnv_segments_{chrom}.png"),
        )

    non_neutral = summarise_cnv(segs_df)
    if not non_neutral.empty:
        print(f"\n  Non-neutral CNV segments:")
        print(non_neutral[["chrom", "start", "end", "mean_log2",
                            "cn_state", "size_mb"]].to_string(index=False))

    return segs_df


if __name__ == "__main__":
    run_all()

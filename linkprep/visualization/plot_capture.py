"""
Visualizations for Capture Hi-C loop calling results.

Produces:
  capture_arc_<chrom>.png      — arc plot: bait → OEF interactions
  capture_score_dist.png       — score distribution (CHiCAGO tiers)
  capture_distance_decay.png   — contact frequency vs genomic distance
  capture_summary.png          — per-bait interaction count summary
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
from analysis.capture_loops import (
    generate_sample_capture_data,
    call_loops,
    filter_loops,
)

# CHiCAGO tier colours (matches browser display convention)
TIER_COLORS = {
    "significant": "#E53935",   # red   score ≥ 5
    "moderate":    "#1E88E5",   # blue  score 3–5
    "weak":        "#90A4AE",   # grey  score < 3 (rarely plotted)
}


# ---------------------------------------------------------------------------
# Arc plot
# ---------------------------------------------------------------------------

def plot_arc(loops_df: pd.DataFrame, baitmap_df: pd.DataFrame,
             chrom: str, out_path: str,
             region_start: int | None = None,
             region_end: int | None = None):
    """
    Draw a 1D arc plot of Capture Hi-C interactions for one chromosome.

    Arcs connect bait midpoints to OEF midpoints; colour and height
    encode CHiCAGO score tier.
    """
    sub = loops_df[loops_df["bait_chr"] == chrom].copy()
    bsub = baitmap_df[baitmap_df["chrom"] == chrom]
    if sub.empty:
        print(f"  No loops on {chrom} — skipping arc plot.")
        return

    if region_start is not None:
        sub = sub[(sub["bait_start"] >= region_start) |
                  (sub["oe_start"] >= region_start)]
        sub = sub[(sub["bait_end"]   <= region_end) |
                  (sub["oe_end"]     <= region_end)]

    fig, ax = plt.subplots(figsize=(14, 5))

    bait_mids = ((bsub["start"] + bsub["end"]) / 2).values

    for _, row in sub.iterrows():
        bait_mid = (row["bait_start"] + row["bait_end"]) / 2
        oe_mid   = (row["oe_start"]   + row["oe_end"]) / 2
        x_center = (bait_mid + oe_mid) / 2
        width    = abs(oe_mid - bait_mid)
        height   = width * 0.4   # arc height proportional to span
        color    = TIER_COLORS.get(str(row["significance"]), "#90A4AE")
        alpha    = min(0.9, 0.3 + row["score"] / 20)
        lw       = 0.8 + row["score"] / 10

        arc = mpatches.Arc(
            (x_center, 0), width, height * 2,
            angle=0, theta1=0, theta2=180,
            color=color, lw=lw, alpha=alpha,
        )
        ax.add_patch(arc)

    # Bait ticks
    for bm in bait_mids:
        ax.axvline(bm / 1e6, color="#37474F", lw=1.2, alpha=0.6, zorder=5)

    ax.set_xlim(
        (sub[["bait_start", "oe_start"]].min().min() - 500_000) / 1e6,
        (sub[["bait_end",   "oe_end"  ]].max().max() + 500_000) / 1e6,
    )
    ax.set_ylim(-0.05, 1)
    ax.set_xlabel(f"Position on {chrom} (Mb)", fontsize=11)
    ax.set_ylabel("")
    ax.set_yticks([])
    ax.set_title(
        f"Capture Hi-C Loops — {chrom}  "
        f"({len(sub)} interactions from {len(bsub)} baits)",
        fontsize=12,
    )

    legend_handles = [
        mpatches.Patch(color=TIER_COLORS["significant"], label="Significant (score ≥ 5)"),
        mpatches.Patch(color=TIER_COLORS["moderate"],    label="Moderate (3–5)"),
    ]
    ax.legend(handles=legend_handles, loc="upper right", fontsize=9)
    ax.spines[["top", "right", "left"]].set_visible(False)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Score distribution
# ---------------------------------------------------------------------------

def plot_score_distribution(loops_df: pd.DataFrame, out_path: str,
                             cutoff: float = 5.0):
    fig, ax = plt.subplots(figsize=(8, 4))

    scores = loops_df["score"].values
    bins   = np.linspace(0, max(scores.max(), cutoff + 5), 50)

    # Colour bars by tier
    for i in range(len(bins) - 1):
        lo, hi = bins[i], bins[i + 1]
        mask   = (scores >= lo) & (scores < hi)
        count  = mask.sum()
        if lo >= cutoff:
            color = TIER_COLORS["significant"]
        elif lo >= 3:
            color = TIER_COLORS["moderate"]
        else:
            color = TIER_COLORS["weak"]
        ax.bar(lo, count, width=hi - lo, color=color,
               edgecolor="white", linewidth=0.3, align="edge")

    ax.axvline(cutoff, color="#37474F", lw=2, ls="--",
               label=f"CHiCAGO cutoff = {cutoff}")
    ax.axvline(3, color="#78909C", lw=1, ls=":",
               label="Moderate cutoff = 3")
    ax.set_xlabel("CHiCAGO Score (−log₁₀ p-value)", fontsize=11)
    ax.set_ylabel("Number of interactions", fontsize=11)
    ax.set_title("Capture Hi-C — Interaction Score Distribution", fontsize=12)
    ax.legend(fontsize=9)

    n_sig = (scores >= cutoff).sum()
    n_mod = ((scores >= 3) & (scores < cutoff)).sum()
    ax.text(0.97, 0.95,
            f"Significant: {n_sig}\nModerate: {n_mod}",
            transform=ax.transAxes, ha="right", va="top",
            fontsize=9, bbox=dict(boxstyle="round", fc="white", alpha=0.8))

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Distance decay
# ---------------------------------------------------------------------------

def plot_distance_decay(loops_df: pd.DataFrame, out_path: str):
    cis = loops_df[loops_df["bait_chr"] == loops_df["oe_chr"]].copy()
    if cis.empty:
        return

    dist_mb = cis["distance_bp"].values / 1e6
    scores  = cis["score"].values

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # Panel 1: mean reads vs distance
    bins_mb = np.linspace(0, dist_mb.max(), 30)
    means, edges, _ = axes[0].hist(
        dist_mb, bins=bins_mb,
        weights=cis["n_reads"].values,
        color="#1E88E5", alpha=0.7, edgecolor="white",
    )
    axes[0].set_xlabel("Genomic distance (Mb)")
    axes[0].set_ylabel("Total reads")
    axes[0].set_title("Read depth vs distance")

    # Panel 2: score vs distance scatter
    sig_mask = cis["significance"] == "significant"
    axes[1].scatter(dist_mb[~sig_mask], scores[~sig_mask],
                    s=8, alpha=0.3, color="#90A4AE", label="Moderate/weak")
    axes[1].scatter(dist_mb[sig_mask], scores[sig_mask],
                    s=12, alpha=0.7, color="#E53935", label="Significant")
    axes[1].axhline(5, color="#37474F", lw=1.5, ls="--", label="Cutoff = 5")
    axes[1].set_xlabel("Genomic distance (Mb)")
    axes[1].set_ylabel("CHiCAGO score")
    axes[1].set_title("Score vs distance")
    axes[1].legend(fontsize=8)

    plt.suptitle("Capture Hi-C — Distance Decay Analysis", fontsize=12)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Per-bait summary
# ---------------------------------------------------------------------------

def plot_per_bait_summary(loops_df: pd.DataFrame, out_path: str,
                           cutoff: float = 5.0):
    sig = loops_df[loops_df["score"] >= cutoff]
    counts = sig.groupby("bait_name").size().sort_values(ascending=False)

    if counts.empty:
        print("  No significant interactions — skipping per-bait summary.")
        return

    fig, ax = plt.subplots(figsize=(max(8, len(counts) * 0.5 + 2), 4))
    colors = plt.cm.RdYlGn(np.linspace(0.3, 0.9, len(counts)))
    bars = ax.bar(range(len(counts)), counts.values,
                  color=colors, edgecolor="white", linewidth=0.4)
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(counts.index, rotation=45, ha="right", fontsize=8)
    ax.set_ylabel("Significant interactions (score ≥ 5)")
    ax.set_title("Capture Hi-C — Interactions per Bait", fontsize=12)
    for bar, val in zip(bars, counts.values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.1,
                str(val), ha="center", va="bottom", fontsize=7)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all():
    os.makedirs(config.FIGURES_DIR, exist_ok=True)

    print("Generating synthetic capture Hi-C data …")
    baitmap_df, rmap_df, contacts_df = generate_sample_capture_data(
        n_baits=20, n_oe_per_bait=50,
        bin_size=config.BIN_SIZE,
    )

    print(f"  {len(baitmap_df)} baits, {len(contacts_df)} raw contacts")
    loops_df  = call_loops(contacts_df, cutoff=5.0,
                            n_total_pairs=500_000)
    loops_df  = filter_loops(loops_df)
    n_sig     = (loops_df["score"] >= 5).sum()
    print(f"  After filtering: {len(loops_df)} interactions, {n_sig} significant (score ≥ 5)")

    for chrom in baitmap_df["chrom"].unique():
        plot_arc(
            loops_df, baitmap_df, chrom,
            os.path.join(config.FIGURES_DIR, f"capture_arc_{chrom}.png"),
        )

    plot_score_distribution(
        loops_df,
        os.path.join(config.FIGURES_DIR, "capture_score_dist.png"),
    )
    plot_distance_decay(
        loops_df,
        os.path.join(config.FIGURES_DIR, "capture_distance_decay.png"),
    )
    plot_per_bait_summary(
        loops_df,
        os.path.join(config.FIGURES_DIR, "capture_summary.png"),
    )

    return loops_df, baitmap_df


if __name__ == "__main__":
    run_all()

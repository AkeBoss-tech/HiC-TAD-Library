"""
Visualize structural variants detected from LinkPrep pairs.

Produces:
  sv_translocation_heatmap.png  — inter-chromosomal pair density
  sv_summary_bar.png            — SV count by type
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
from analysis.sv_detection import detect_svs


# ---------------------------------------------------------------------------
# Translocation heatmap
# ---------------------------------------------------------------------------

def plot_translocation_heatmap(translocations: pd.DataFrame, out_path: str,
                                bin_size: int = 100_000):
    if translocations.empty:
        print("  No translocations to plot.")
        return

    chroms = sorted(set(translocations["chrom1"]) | set(translocations["chrom2"]))
    n_chroms = len(chroms)
    chrom_idx = {c: i for i, c in enumerate(chroms)}

    fig, ax = plt.subplots(figsize=(8, 6))

    for _, row in translocations.iterrows():
        i = chrom_idx[row["chrom1"]]
        j = chrom_idx[row["chrom2"]]
        ax.scatter(j, i, s=row["n_pairs"] * 5, alpha=0.6,
                   color="#E91E63", edgecolors="white", linewidths=0.3)

    ax.set_xticks(range(n_chroms))
    ax.set_xticklabels(chroms, rotation=45, ha="right")
    ax.set_yticks(range(n_chroms))
    ax.set_yticklabels(chroms)
    ax.set_title("LinkPrep — Translocation Candidates\n"
                 "(bubble size ∝ supporting pairs)", fontsize=11)
    ax.set_xlabel("Chromosome 2")
    ax.set_ylabel("Chromosome 1")
    ax.grid(True, alpha=0.2)

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# SV summary bar chart
# ---------------------------------------------------------------------------

def plot_sv_summary(svs: dict, out_path: str):
    types  = ["translocations", "inversions", "deletions"]
    labels = ["Translocations", "Inversions", "Deletions"]
    counts = [len(svs.get(t, pd.DataFrame())) for t in types]
    colors = ["#E91E63", "#FF9800", "#9C27B0"]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, counts, color=colors, edgecolor="white")
    for bar, cnt in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + 0.3,
                str(cnt), ha="center", va="bottom", fontsize=11, fontweight="bold")
    ax.set_ylabel("Candidate bins")
    ax.set_title("LinkPrep — SV Detection Summary")
    ax.set_ylim(0, max(counts) * 1.2 + 1)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Inversion lollipop
# ---------------------------------------------------------------------------

def plot_inversions(inversions: pd.DataFrame, out_path: str):
    if inversions.empty:
        print("  No inversions to plot.")
        return

    for chrom, grp in inversions.groupby("chrom"):
        fig, ax = plt.subplots(figsize=(10, 3))
        for _, row in grp.iterrows():
            x1 = row["bin1"] / 1e6
            x2 = row["bin2"] / 1e6
            ax.plot([x1, x2], [0, 0], lw=2 + row["n_pairs"] * 0.2,
                    color="#FF9800", solid_capstyle="round")
            ax.scatter([x1, x2], [0, 0], s=30, color="#E65100", zorder=5)

        ax.set_xlabel(f"Position on {chrom} (Mb)")
        ax.set_yticks([])
        ax.set_title(f"Putative inversions — {chrom}")
        plt.tight_layout()
        fname = out_path.replace(".png", f"_{chrom}.png")
        plt.savefig(fname, dpi=150)
        plt.close()
        print(f"  Saved → {fname}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all(pairs_path: str):
    os.makedirs(config.FIGURES_DIR, exist_ok=True)
    svs = detect_svs(pairs_path, bin_size=100_000, min_pairs=5)

    plot_translocation_heatmap(
        svs["translocations"],
        os.path.join(config.FIGURES_DIR, "sv_translocation_heatmap.png"),
    )
    plot_sv_summary(
        svs,
        os.path.join(config.FIGURES_DIR, "sv_summary.png"),
    )
    plot_inversions(
        svs["inversions"],
        os.path.join(config.FIGURES_DIR, "sv_inversions.png"),
    )
    return svs


if __name__ == "__main__":
    if not os.path.exists(config.PAIRS_FILE):
        sys.exit(f"Pairs file not found: {config.PAIRS_FILE}\n"
                 "Run: python sample_data/generate_sample_data.py")
    run_all(config.PAIRS_FILE)

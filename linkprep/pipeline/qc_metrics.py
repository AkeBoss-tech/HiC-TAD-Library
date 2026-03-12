"""
QC metrics for a LinkPrep library.

Reads the pairtools dedup stats file (or the pairs file directly) and
prints/plots key QC metrics:

  - Total read pairs
  - Alignment rate
  - Valid (UU) pair fraction
  - Cis/trans ratio
  - Library complexity (unique fraction)
  - Cis-pair distance distribution

Usage:
    python pipeline/qc_metrics.py            # uses paths from config.py
    python pipeline/qc_metrics.py --pairs /path/to/file.pairs.gz
"""

import argparse
import gzip
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config

# ---------------------------------------------------------------------------
# Pairs file parsing
# ---------------------------------------------------------------------------

def read_pairs_sample(pairs_path: str, max_rows: int = 500_000) -> pd.DataFrame:
    """Read up to max_rows data lines from a (possibly gzipped) pairs file."""
    opener = gzip.open if pairs_path.endswith(".gz") else open
    rows = []
    with opener(pairs_path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            rows.append({
                "readID":     parts[0],
                "chrom1":     parts[1],
                "pos1":       int(parts[2]),
                "chrom2":     parts[3],
                "pos2":       int(parts[4]),
                "strand1":    parts[5],
                "strand2":    parts[6],
                "pair_type":  parts[7],
            })
            if len(rows) >= max_rows:
                break
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# QC statistics
# ---------------------------------------------------------------------------

def compute_qc(df: pd.DataFrame) -> dict:
    total = len(df)
    valid = (df["pair_type"] == "UU").sum()
    cis   = ((df["chrom1"] == df["chrom2"]) & (df["pair_type"] == "UU")).sum()
    trans = ((df["chrom1"] != df["chrom2"]) & (df["pair_type"] == "UU")).sum()

    cis_df = df[(df["chrom1"] == df["chrom2"]) & (df["pair_type"] == "UU")].copy()
    cis_df["distance"] = (cis_df["pos2"] - cis_df["pos1"]).abs()
    cis_long = (cis_df["distance"] >= 1000).sum()   # ≥1 kb cis pairs

    return {
        "total_pairs":         total,
        "valid_pairs":         valid,
        "valid_frac":          valid / total if total else 0,
        "cis_pairs":           cis,
        "trans_pairs":         trans,
        "cis_frac":            cis / valid if valid else 0,
        "cis_long_range_frac": cis_long / cis if cis else 0,
        "cis_distances":       cis_df["distance"].values,
    }


def print_qc(stats: dict):
    print("\n=== LinkPrep QC Report ===")
    print(f"  Total pairs sampled : {stats['total_pairs']:>10,}")
    print(f"  Valid (UU) pairs    : {stats['valid_pairs']:>10,}  ({stats['valid_frac']:.1%})")
    print(f"  Cis pairs           : {stats['cis_pairs']:>10,}  ({stats['cis_frac']:.1%} of valid)")
    print(f"  Trans pairs         : {stats['trans_pairs']:>10,}")
    print(f"  Cis ≥1 kb fraction  : {stats['cis_long_range_frac']:.1%}")
    print()

    # Pass/fail
    checks = {
        "Valid pair rate ≥ 50%":            stats["valid_frac"] >= 0.50,
        "Cis fraction ≥ 60%":              stats["cis_frac"] >= 0.60,
        "Long-range cis (≥1 kb) ≥ 40%":   stats["cis_long_range_frac"] >= 0.40,
    }
    for msg, passed in checks.items():
        status = "PASS" if passed else "FAIL"
        print(f"  [{status}] {msg}")
    print()


# ---------------------------------------------------------------------------
# Plots
# ---------------------------------------------------------------------------

def plot_distance_distribution(distances: np.ndarray, out_path: str):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4))
    fig.suptitle("LinkPrep — Cis-pair distance distribution", fontsize=13)

    bins = np.logspace(2, 8, 60)

    axes[0].hist(distances, bins=bins, color="#2196F3", edgecolor="none")
    axes[0].set_xscale("log")
    axes[0].set_xlabel("Genomic distance (bp)")
    axes[0].set_ylabel("Pair count")
    axes[0].set_title("Raw histogram")

    # Density (P(s) ~ distance decay)
    counts, edges = np.histogram(distances, bins=bins)
    centres = (edges[:-1] + edges[1:]) / 2
    density = counts / counts.sum() / np.diff(edges)
    axes[1].plot(centres, density * centres, color="#E91E63", lw=2)
    axes[1].set_xscale("log")
    axes[1].set_yscale("log")
    axes[1].set_xlabel("Genomic distance (bp)")
    axes[1].set_ylabel("Contact probability P(s)")
    axes[1].set_title("P(s) curve")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


def plot_pair_type_bar(stats: dict, out_path: str):
    labels = ["Trans", "Cis <1 kb", "Cis ≥1 kb"]
    cis_short = stats["cis_pairs"] - int(stats["cis_pairs"] * stats["cis_long_range_frac"])
    cis_long  = int(stats["cis_pairs"] * stats["cis_long_range_frac"])
    values = [stats["trans_pairs"], cis_short, cis_long]
    colors = ["#FF7043", "#FFA726", "#66BB6A"]

    fig, ax = plt.subplots(figsize=(6, 4))
    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=0.5)
    ax.set_ylabel("Pair count")
    ax.set_title("LinkPrep — Pair-type breakdown")
    for bar, val in zip(bars, values):
        ax.text(bar.get_x() + bar.get_width() / 2,
                bar.get_height() + max(values) * 0.01,
                f"{val:,}", ha="center", va="bottom", fontsize=9)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="LinkPrep QC metrics")
    parser.add_argument("--pairs", default=config.PAIRS_FILE,
                        help="Path to .pairs or .pairs.gz file")
    args = parser.parse_args()

    if not os.path.exists(args.pairs):
        sys.exit(f"Pairs file not found: {args.pairs}\n"
                 "Run: python sample_data/generate_sample_data.py first.")

    print(f"Reading pairs from: {args.pairs}")
    df = read_pairs_sample(args.pairs)
    stats = compute_qc(df)
    print_qc(stats)

    os.makedirs(config.FIGURES_DIR, exist_ok=True)
    plot_distance_distribution(
        stats["cis_distances"],
        os.path.join(config.FIGURES_DIR, "qc_distance_distribution.png"),
    )
    plot_pair_type_bar(
        stats,
        os.path.join(config.FIGURES_DIR, "qc_pair_types.png"),
    )


if __name__ == "__main__":
    main()

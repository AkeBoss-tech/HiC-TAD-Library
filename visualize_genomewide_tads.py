#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def _maybe_numeric_columns(df: pd.DataFrame, suffix: str) -> list[str]:
    cols = [c for c in df.columns if c.endswith(suffix)]
    return [c for c in cols if pd.api.types.is_numeric_dtype(df[c])]


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Visualize genome-wide TAD summary tables.")
    p.add_argument("--prefix", required=True, help="Prefix stem, e.g. data/processed/genomewide_tads/mouse_microc_100kb")
    return p


def main() -> None:
    args = build_parser().parse_args()
    prefix = Path(args.prefix)

    tads = pd.read_csv(prefix.with_suffix(".tads.tsv"), sep="\t")
    tad_summary = pd.read_csv(prefix.with_suffix(".tad_summary.tsv"), sep="\t")
    boundary_summary = pd.read_csv(prefix.with_suffix(".boundary_summary.tsv"), sep="\t")

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    axes[0, 0].hist(tads["length_bp"] / 1e6, bins=40, color="#b44f3a", edgecolor="white")
    axes[0, 0].set_title("TAD Size Distribution")
    axes[0, 0].set_xlabel("TAD size (Mb)")
    axes[0, 0].set_ylabel("Count")

    chrom_counts = tads["chrom"].value_counts().sort_index()
    axes[0, 1].bar(chrom_counts.index, chrom_counts.values, color="#3e6c8f")
    axes[0, 1].set_title("TADs per Chromosome")
    axes[0, 1].set_ylabel("Count")
    axes[0, 1].tick_params(axis="x", rotation=90)

    assay_cols = _maybe_numeric_columns(tad_summary, "_mean")
    assay_cols = [c for c in assay_cols if c.split("_")[0] in {"atac", "chip", "rna"}]
    if assay_cols:
        plot_df = tad_summary[assay_cols].copy()
        plot_df.columns = [c.replace("_mean", "") for c in plot_df.columns]
        axes[1, 0].boxplot([plot_df[c] for c in plot_df.columns], labels=plot_df.columns)
        axes[1, 0].set_title("Assay Signal per TAD")
        axes[1, 0].set_ylabel("Overlap-weighted mean signal")
    else:
        axes[1, 0].text(0.5, 0.5, "No assay summary columns found", ha="center", va="center")
        axes[1, 0].set_title("Assay Signal per TAD")
        axes[1, 0].set_axis_off()

    boundary_cols = _maybe_numeric_columns(boundary_summary, "_mean")
    boundary_cols = [c for c in boundary_cols if c.split("_")[0] in {"atac", "chip", "rna"}]
    if "chip_boundary_mean" in boundary_summary.columns and "rna_boundary_mean" in boundary_summary.columns:
        axes[1, 1].scatter(
            boundary_summary["chip_boundary_mean"],
            boundary_summary["rna_boundary_mean"],
            s=8,
            alpha=0.4,
            color="#2c7a7b",
        )
        axes[1, 1].set_xlabel("Boundary ChIP mean")
        axes[1, 1].set_ylabel("Boundary RNA mean")
        axes[1, 1].set_title("Boundary ChIP vs RNA")
    elif boundary_cols:
        first = boundary_cols[0]
        axes[1, 1].hist(boundary_summary[first], bins=40, color="#6c8f3e", edgecolor="white")
        axes[1, 1].set_title(f"{first.replace('_mean', '')} boundary signal")
        axes[1, 1].set_xlabel("Signal")
    else:
        axes[1, 1].text(0.5, 0.5, "No boundary assay summary columns found", ha="center", va="center")
        axes[1, 1].set_title("Boundary Assay Summary")
        axes[1, 1].set_axis_off()

    fig.suptitle(prefix.name.replace("_", " "))
    fig.tight_layout()

    out = prefix.with_suffix(".summary.png")
    fig.savefig(out, dpi=180, bbox_inches="tight")
    print(f"Saved {out}")


if __name__ == "__main__":
    main()

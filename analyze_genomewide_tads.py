#!/usr/bin/env python3
from __future__ import annotations

import argparse
from pathlib import Path

from src.genomewide_tads import (
    annotate_tads_with_assays,
    call_genomewide_tads,
    load_interval_signal_table,
)
from src.loaders import load_cooler


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Call genome-wide TADs from a cooler contact map and annotate them with ATAC/ChIP/RNA interval tracks.",
    )
    parser.add_argument("--cooler", required=True, help="Filename in data/raw, for example mouse_microc.mcool")
    parser.add_argument("--resolution", type=int, default=5_000, help="Cooler resolution in base pairs")
    parser.add_argument("--window-bp", type=int, default=100_000, help="Insulation window for genome-wide TAD calling")
    parser.add_argument("--ignore-diags", type=int, default=2, help="Diagonals to ignore in insulation")
    parser.add_argument("--boundary-class-min", default="weak", choices=["sub_threshold", "weak", "strong"])
    parser.add_argument("--min-chrom-size-bp", type=int, default=500_000)
    parser.add_argument("--min-tad-length-bp", type=int, default=100_000)
    parser.add_argument("--max-tad-length-bp", type=int, default=3_000_000)
    parser.add_argument("--chromosomes", nargs="*", default=None, help="Optional chromosome allow-list")
    parser.add_argument("--atac", help="BED/BEDGRAPH/TSV interval file with ATAC signal")
    parser.add_argument("--chip", help="BED/BEDGRAPH/TSV interval file with ChIP signal, e.g. CTCF")
    parser.add_argument("--rna", help="BED/BEDGRAPH/TSV interval file with RNA signal over gene intervals")
    parser.add_argument("--boundary-window-bp", type=int, default=25_000)
    parser.add_argument("--outdir", default="data/processed/genomewide_tads")
    parser.add_argument("--prefix", default="genomewide")
    return parser


def main() -> None:
    args = build_parser().parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    clr = load_cooler(args.cooler, resolution=args.resolution)
    insulation_df, boundary_df, tad_df = call_genomewide_tads(
        clr,
        window_bp=args.window_bp,
        chromosomes=args.chromosomes,
        ignore_diags=args.ignore_diags,
        boundary_class_min=args.boundary_class_min,
        min_chrom_size_bp=args.min_chrom_size_bp,
        min_tad_length_bp=args.min_tad_length_bp,
        max_tad_length_bp=args.max_tad_length_bp,
    )

    atac_df = load_interval_signal_table(args.atac) if args.atac else None
    chip_df = load_interval_signal_table(args.chip) if args.chip else None
    rna_df = load_interval_signal_table(args.rna) if args.rna else None

    tad_summary_df, boundary_summary_df = annotate_tads_with_assays(
        tad_df,
        chromsizes=clr.chromsizes,
        atac_df=atac_df,
        chip_df=chip_df,
        rna_df=rna_df,
        boundary_window_bp=args.boundary_window_bp,
    )

    stem = outdir / args.prefix
    insulation_df.to_csv(f"{stem}.insulation.tsv", sep="\t", index=False)
    boundary_df.to_csv(f"{stem}.boundaries.tsv", sep="\t", index=False)
    tad_df.to_csv(f"{stem}.tads.tsv", sep="\t", index=False)
    tad_summary_df.to_csv(f"{stem}.tad_summary.tsv", sep="\t", index=False)
    boundary_summary_df.to_csv(f"{stem}.boundary_summary.tsv", sep="\t", index=False)

    print(f"Saved insulation track to {stem}.insulation.tsv")
    print(f"Saved boundary table to {stem}.boundaries.tsv")
    print(f"Saved TAD calls to {stem}.tads.tsv")
    print(f"Saved TAD assay summary to {stem}.tad_summary.tsv")
    print(f"Saved boundary assay summary to {stem}.boundary_summary.tsv")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Regulatory Activity Pipeline.

Examples:
    python pipeline/stage_regulatory_activity.py \
        --organism mouse --assembly mm10 --cell-type CL:0000540 \
        --interval chr12:26000000-28000000
"""

from __future__ import annotations

import argparse
import json
import os
import sys
from pathlib import Path

# Add repo root to sys.path so src.* imports work
REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if REPO_ROOT not in sys.path:
    sys.path.insert(0, REPO_ROOT)

from src.regulatory_activity import (  # noqa: E402
    AlphaGenomePredictionAdapter,
    LocusRequest,
    build_locus_context,
    build_manifest,
    compare_reference_vs_edit,
    load_external_evidence,
    rank_regulatory_activity,
    save_outputs,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Rank regulatory activity from DNA sequence + 3D context.")
    parser.add_argument("--organism", required=True, choices=["mouse", "human"])
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--cell-type", required=True)
    parser.add_argument("--interval", help="Genomic interval, e.g. chr12:26000000-28000000")
    parser.add_argument("--anchor-interval", help="Required with --sequence-fasta when the sequence is synthetic.")
    parser.add_argument("--sequence-fasta", help="Optional FASTA file containing a synthetic or edited sequence.")
    parser.add_argument("--compare-edit", help="Optional JSON file describing deletion/substitution/insertion edit.")
    parser.add_argument("--fasta-path", help="Optional reference FASTA override.")
    parser.add_argument("--annotation-path", help="Optional local annotation file override (reserved for future use).")
    parser.add_argument("--cooler-filename", help="Optional local .mcool filename under data/raw for observed contact support.")
    parser.add_argument("--cooler-resolution", type=int, default=5000)
    parser.add_argument("--window-bp", type=int, default=1_048_576)
    parser.add_argument("--ranking-resolution-bp", type=int, default=2_048)
    return parser.parse_args()


def _read_fasta_sequence(path: str) -> str:
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    return "".join(lines).upper()


def main(mode: str = "single") -> dict:
    args = parse_args()
    sequence = _read_fasta_sequence(args.sequence_fasta) if args.sequence_fasta else None
    compare_edit = None
    if args.compare_edit:
        with open(args.compare_edit) as f:
            compare_edit = json.load(f)

    request = LocusRequest(
        organism=args.organism,
        assembly=args.assembly,
        cell_type=args.cell_type,
        interval=args.interval,
        sequence=sequence,
        anchor_interval=args.anchor_interval,
        compare_edit=compare_edit,
        window_bp=args.window_bp,
        ranking_resolution_bp=args.ranking_resolution_bp,
        fasta_path=args.fasta_path,
        annotation_path=args.annotation_path,
        local_cooler_filename=args.cooler_filename,
        local_cooler_resolution=args.cooler_resolution,
    )

    context = build_locus_context(request)
    adapter = AlphaGenomePredictionAdapter()
    prediction = adapter.predict_reference(context)
    external_evidence, external_provenance = load_external_evidence(context)
    reference_bins, promoter_links, gene_summaries, scoring_summary = rank_regulatory_activity(
        context,
        prediction,
        external_evidence,
    )

    edited_bins = None
    edited_links = None
    edited_gene_summaries = None
    edited_prediction = None
    perturbation_summary = None

    if request.compare_edit:
        edited_bins, edited_links, edited_gene_summaries, _, perturbation_summary, edited_prediction = compare_reference_vs_edit(
            context,
            adapter,
            reference_bins,
            request.compare_edit,
            external_evidence,
        )

    manifest = build_manifest(
        context,
        reference_bins,
        promoter_links,
        gene_summaries,
        prediction,
        external_provenance,
        scoring_summary,
        edited_bins=edited_bins,
        edited_links=edited_links,
        edited_gene_summaries=edited_gene_summaries,
        edited_prediction=edited_prediction,
        perturbation_summary=perturbation_summary,
    )
    outputs = save_outputs(context, manifest, reference_bins, gene_summaries, edited_bins=edited_bins)
    print("Regulatory activity analysis complete.")
    for label, path in outputs.items():
        print(f"  {label}: {path}")
    return manifest


if __name__ == "__main__":
    main()

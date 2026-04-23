from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any

from src.hg_dt.models.alphagenome import CELL_TYPES
from src.regulatory_activity import (
    AlphaGenomePredictionAdapter,
    LocusRequest,
    build_locus_context,
    build_manifest,
    compare_reference_vs_edit,
    load_external_evidence,
    manifest_slug,
    rank_regulatory_activity,
    save_outputs,
)


GENE_PRESETS = {
    "Sox11": {"organism": "mouse", "assembly": "mm10", "interval": "chr12:26000000-28000000", "cell_type": "CL:0000540"},
    "Mir9-2": {"organism": "mouse", "assembly": "mm10", "interval": "chr13:83500000-84500000", "cell_type": "CL:0000540"},
    "TAL1": {"organism": "human", "assembly": "hg38", "interval": "chr1:47210000-47260000", "cell_type": "Jurkat"},
    "OCT4": {"organism": "human", "assembly": "hg38", "interval": "chr6:30630000-31630000", "cell_type": "H1"},
}


def project_root() -> Path:
    return Path(__file__).resolve().parent.parent


def media_root() -> Path:
    return project_root() / "media"


def processed_root() -> Path:
    return project_root() / "data" / "processed"


def list_saved_runs() -> list[dict[str, Any]]:
    runs = []
    for path in sorted(processed_root().glob("regulatory_activity_*.json"), reverse=True):
        try:
            with open(path) as f:
                payload = json.load(f)
            slug = path.stem.replace("regulatory_activity_", "", 1)
            runs.append(
                {
                    "slug": slug,
                    "path": str(path),
                    "organism": payload.get("input", {}).get("organism"),
                    "assembly": payload.get("input", {}).get("assembly"),
                    "cell_type": payload.get("input", {}).get("cell_type"),
                    "interval": payload.get("input", {}).get("interval") or payload.get("input", {}).get("anchor_interval"),
                    "window": payload.get("input", {}).get("window"),
                    "candidate_count": len(payload.get("candidate_elements", [])),
                    "gene_count": len(payload.get("gene_summaries", [])),
                }
            )
        except Exception:
            continue
    return runs


def load_run(slug: str) -> dict[str, Any]:
    path = processed_root() / f"regulatory_activity_{slug}.json"
    if not path.exists():
        raise FileNotFoundError(f"No run found for slug {slug}")
    with open(path) as f:
        return json.load(f)


def available_cell_types() -> list[str]:
    return sorted(CELL_TYPES.keys())


def _read_fasta_sequence(path: str) -> str:
    lines = []
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                continue
            lines.append(line.strip())
    return "".join(lines).upper()


def _maybe_edit_payload(form: dict[str, Any]) -> dict[str, Any] | None:
    mode = str(form.get("edit_mode") or "").strip().lower()
    if not mode or mode == "none":
        return None
    start = form.get("edit_start")
    end = form.get("edit_end")
    if not start:
        return None
    payload = {
        "mode": mode,
        "start": int(start),
        "end": int(end) if end else int(start),
    }
    seq = str(form.get("edit_sequence") or "").strip().upper()
    if seq:
        payload["sequence"] = seq
    return payload


def readiness_snapshot() -> dict[str, Any]:
    return {
        "alpha_key": bool(os.environ.get("ALPHA_GENOME_API_KEY")),
        "nvidia_key": bool(os.environ.get("NVIDIA_API_KEY") or os.environ.get("NVIDIA_ESM_FOLD_API_KEY")),
        "mouse_fasta": (project_root() / "data" / "raw" / "mm10.fa").exists(),
        "mouse_microc": (project_root() / "data" / "raw" / "mouse_microc.mcool").exists(),
        "human_fasta": (project_root() / "data" / "references" / "hg38.fa").exists(),
        "gencode_gtf": (project_root() / "data" / "references" / "gencode.v49.annotation.gtf").exists(),
    }


def run_regulatory_activity(form: dict[str, Any]) -> dict[str, Any]:
    organism = str(form["organism"]).strip().lower()
    assembly = str(form["assembly"]).strip()
    cell_type = str(form["cell_type"]).strip()
    interval = str(form.get("interval") or "").strip() or None
    anchor_interval = str(form.get("anchor_interval") or "").strip() or None
    sequence_fasta = str(form.get("sequence_fasta") or "").strip() or None
    sequence = _read_fasta_sequence(sequence_fasta) if sequence_fasta else None
    compare_edit = _maybe_edit_payload(form)
    cooler_filename = str(form.get("cooler_filename") or "").strip() or None
    cooler_resolution = int(form.get("cooler_resolution") or 5000)
    ranking_resolution = int(form.get("ranking_resolution_bp") or 2048)
    window_bp = int(form.get("window_bp") or 1_048_576)

    request = LocusRequest(
        organism=organism,
        assembly=assembly,
        cell_type=cell_type,
        interval=interval,
        anchor_interval=anchor_interval,
        sequence=sequence,
        compare_edit=compare_edit,
        window_bp=window_bp,
        ranking_resolution_bp=ranking_resolution,
        fasta_path=str(form.get("fasta_path") or "").strip() or None,
        local_cooler_filename=cooler_filename,
        local_cooler_resolution=cooler_resolution,
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

    if compare_edit:
        edited_bins, edited_links, edited_gene_summaries, _, perturbation_summary, edited_prediction = compare_reference_vs_edit(
            context,
            adapter,
            reference_bins,
            compare_edit,
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
    slug = manifest_slug(context)
    return {
        "slug": slug,
        "manifest": manifest,
        "outputs": outputs,
    }

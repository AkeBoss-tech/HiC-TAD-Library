#!/usr/bin/env python3
"""
Stage 0 — Variant-context bridge between synthesis and protein pipeline.

Builds a lightweight JSON manifest that carries forward the strongest synthesis
signal (boundary collapse at the Unc5b insulator) and, when possible, enriches
that context with Ensembl transcript/translation metadata for differential runs.
"""

from __future__ import annotations

import json
import os
import urllib.error
import urllib.request
from datetime import datetime, timezone

import numpy as np

GENE_SYMBOL = "Unc5b"
SPECIES = "mus_musculus"
WT_CACHE = "data/processed/synthesis_unc5b_wt_cache.npz"
OUT_JSON = "data/processed/pipeline_variant_context.json"

# Repository-grounded prior from synthesis Exp 2 / synthesis_report.py.
BOUNDARY_COLLAPSE_DELTA = 0.268


def load_wt_cache_summary(cache_path: str = WT_CACHE) -> dict:
    if not os.path.exists(cache_path):
        return {"available": False, "path": cache_path}

    data = np.load(cache_path)
    contact = data["contact_mat"]
    ctcf = data["ctcf_signal"]
    return {
        "available": True,
        "path": cache_path,
        "contact_map_shape": list(contact.shape),
        "ctcf_bins": int(len(ctcf)),
        "resolution_bp": float(data["resolution"]),
    }


def fetch_json(url: str) -> object:
    req = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "Content-Type": "application/json"},
        method="GET",
    )
    with urllib.request.urlopen(req, timeout=20) as resp:
        return json.loads(resp.read().decode())


def fetch_ensembl_isoforms(species: str = SPECIES, gene_symbol: str = GENE_SYMBOL) -> tuple[list[dict], str | None]:
    try:
        xref_url = (
            f"https://rest.ensembl.org/xrefs/symbol/{species}/{gene_symbol}"
            "?object_type=gene"
        )
        xrefs = fetch_json(xref_url)
        if not xrefs:
            return [], "No Ensembl gene match found"

        gene_id = xrefs[0]["id"]
        lookup_url = f"https://rest.ensembl.org/lookup/id/{gene_id}?expand=1"
        gene_info = fetch_json(lookup_url)

        isoforms = []
        for transcript in gene_info.get("Transcript", []):
            translation = transcript.get("Translation")
            if not translation or "id" not in translation:
                continue

            protein_id = translation["id"]
            seq_url = f"https://rest.ensembl.org/sequence/id/{protein_id}?type=protein"
            seq_info = fetch_json(seq_url)
            sequence = seq_info.get("seq", "")
            if not sequence:
                continue

            isoforms.append(
                {
                    "transcript_id": transcript.get("id"),
                    "translation_id": protein_id,
                    "display_name": transcript.get("display_name") or transcript.get("id"),
                    "is_canonical": bool(transcript.get("is_canonical")),
                    "protein_length": len(sequence),
                    "protein_sequence": sequence,
                }
            )

        isoforms.sort(key=lambda row: (not row["is_canonical"], -row["protein_length"], row["display_name"]))
        return isoforms, None
    except (urllib.error.URLError, TimeoutError, KeyError, json.JSONDecodeError) as exc:
        return [], f"Ensembl lookup failed: {exc}"


def build_context() -> dict:
    wt_cache = load_wt_cache_summary()
    isoforms, isoform_warning = fetch_ensembl_isoforms()

    context = {
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "gene": {
            "symbol": GENE_SYMBOL,
            "species": SPECIES,
        },
        "wildtype_cache": wt_cache,
        "boundary_collapse": {
            "variant_id": "del_4ctcf",
            "delta_insulation": BOUNDARY_COLLAPSE_DELTA,
            "interpretation": "Major collapse inferred from synthesis Experiment 2",
            "expression_prior": {
                "direction": "up",
                "fold_change": None,
                "source": "inference_from_boundary_collapse",
            },
        },
        "isoforms": isoforms,
        "warnings": [],
        "fallback_variant": {
            "id": "boundary_collapse_proxy",
            "label": "Boundary-collapse proxy truncation",
            "strategy": "truncate_c_terminal_tail",
            "truncate_residues": 16,
            "reason": "Used when no alternate translated isoform sequence is available",
        },
    }

    if isoform_warning:
        context["warnings"].append(isoform_warning)
    if not wt_cache["available"]:
        context["warnings"].append(
            "WT synthesis cache missing; Stage 1 can still run, but chromatin context will be weaker."
        )

    return context


def main() -> dict:
    context = build_context()
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as f:
        json.dump(context, f, indent=2)
    print(f"Variant context saved -> {OUT_JSON}")
    if context["warnings"]:
        print("Warnings:")
        for warning in context["warnings"]:
            print(f"  - {warning}")
    return context


if __name__ == "__main__":
    main()

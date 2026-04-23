from __future__ import annotations

import traceback
from pathlib import Path

from flask import Blueprint, Response, current_app, flash, jsonify, redirect, render_template, request, send_file, url_for

from flask_dashboard.services import (
    GENE_PRESETS,
    available_cell_types,
    list_saved_runs,
    load_run,
    media_root,
    readiness_snapshot,
    run_regulatory_activity,
)

bp = Blueprint("dashboard", __name__)


def _preset_payload(name: str) -> dict:
    payload = GENE_PRESETS.get(name, {}).copy()
    payload.setdefault("window_bp", 1_048_576)
    payload.setdefault("ranking_resolution_bp", 2_048)
    payload.setdefault("cooler_resolution", 5_000)
    payload.setdefault("cooler_filename", "mouse_microc.mcool" if payload.get("organism") == "mouse" else "")
    payload.setdefault("anchor_interval", "")
    payload.setdefault("sequence_fasta", "")
    payload.setdefault("edit_mode", "none")
    payload.setdefault("edit_start", "")
    payload.setdefault("edit_end", "")
    payload.setdefault("edit_sequence", "")
    payload.setdefault("fasta_path", "")
    return payload


@bp.get("/")
def index():
    preset_name = request.args.get("preset", "Sox11")
    defaults = _preset_payload(preset_name)
    return render_template(
        "index.html",
        title="Genome Control Room",
        readiness=readiness_snapshot(),
        presets=GENE_PRESETS,
        selected_preset=preset_name,
        defaults=defaults,
        saved_runs=list_saved_runs()[:6],
        cell_types=available_cell_types(),
    )


@bp.get("/saved")
def saved():
    return render_template(
        "saved_runs.html",
        title="Saved Analyses",
        saved_runs=list_saved_runs(),
    )


@bp.post("/analyze")
def analyze():
    form_data = request.form.to_dict()
    try:
        result = run_regulatory_activity(form_data)
        return redirect(url_for("dashboard.result", slug=result["slug"]))
    except Exception as exc:
        flash(f"Analysis failed: {exc}", "error")
        flash(traceback.format_exc(), "trace")
        preset_name = request.form.get("preset_name", "Sox11")
        defaults = _preset_payload(preset_name)
        defaults.update(form_data)
        return render_template(
            "index.html",
            title="Genome Control Room",
            readiness=readiness_snapshot(),
            presets=GENE_PRESETS,
            selected_preset=preset_name,
            defaults=defaults,
            saved_runs=list_saved_runs()[:6],
            cell_types=available_cell_types(),
        ), 500


@bp.get("/result/<slug>")
def result(slug: str):
    payload = load_run(slug)
    outputs = {
        "overview": f"regulatory_activity_{slug}_overview.png",
        "top_candidates": f"regulatory_activity_{slug}_top_candidates.png",
        "delta": f"regulatory_activity_{slug}_delta.png",
    }
    existing_outputs = {name: path for name, path in outputs.items() if (media_root() / path).exists()}
    top_candidates = payload.get("candidate_elements", [])[:15]
    top_genes = payload.get("gene_summaries", [])[:12]
    return render_template(
        "result.html",
        title=f"Analysis {slug}",
        slug=slug,
        payload=payload,
        outputs=existing_outputs,
        top_candidates=top_candidates,
        top_genes=top_genes,
    )


@bp.get("/api/result/<slug>")
def result_api(slug: str):
    return jsonify(load_run(slug))


@bp.get("/media/<path:filename>")
def media_file(filename: str):
    path = media_root() / filename
    if not path.exists():
        return Response("Not found", status=404)
    return send_file(path)


@bp.get("/download/<slug>")
def download_manifest(slug: str):
    project_root = Path(current_app.config["PROJECT_ROOT"])
    path = project_root / "data" / "processed" / f"regulatory_activity_{slug}.json"
    if not path.exists():
        return Response("Not found", status=404)
    return send_file(path, as_attachment=True, download_name=path.name)

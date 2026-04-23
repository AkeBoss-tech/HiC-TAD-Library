from __future__ import annotations

import json
import math
import os
import textwrap
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from src.hg_dt.models.alphagenome import AlphaGenomeConnector
from src.loaders import load_cooler
from src.tad_boundaries import compute_insulation_delta, parse_coordinates

DEFAULT_WINDOW_BP = 1_048_576
DEFAULT_BIN_BP = 2_048
ENSEMBL_SPECIES = {
    "mouse": "mus_musculus",
    "human": "homo_sapiens",
}


@dataclass
class LocusRequest:
    organism: str
    assembly: str
    cell_type: str
    interval: Optional[str] = None
    sequence: Optional[str] = None
    anchor_interval: Optional[str] = None
    compare_edit: Optional[dict[str, Any]] = None
    window_bp: int = DEFAULT_WINDOW_BP
    ranking_resolution_bp: int = DEFAULT_BIN_BP
    fasta_path: Optional[str] = None
    annotation_path: Optional[str] = None
    local_cooler_filename: Optional[str] = None
    local_cooler_resolution: int = 5_000


@dataclass
class CandidateElement:
    chrom: str
    start: int
    end: int
    bin_index: int
    activity_score: float
    sequence_activity_score: float
    promoter_link_score: float
    contact_score: float
    boundary_context_score: float
    cell_type_support_score: float
    promoter_activity_score: float
    ctcf_score: float
    observed_contact_score: float
    d3_domain_id: int
    element_type: str
    linked_gene_ids: list[str] = field(default_factory=list)
    linked_gene_symbols: list[str] = field(default_factory=list)
    evidence_sources: list[str] = field(default_factory=list)
    delta_activity_score: Optional[float] = None


@dataclass
class GeneActivitySummary:
    gene_id: str
    gene_symbol: str
    chrom: str
    start: int
    end: int
    domain_id: int
    promoter_bin: int
    promoter_activity_score: float
    nearest_elements: list[int]
    supporting_elements: list[int]
    linked_element_bins: list[int]
    predicted_expression_change: float
    predicted_splice_change: float
    linked_element_score: float
    evidence_sources: list[str] = field(default_factory=list)


@dataclass
class LocusContext:
    request: LocusRequest
    chrom: str
    interval_start: int
    interval_end: int
    window_start: int
    window_end: int
    query_start: int
    query_end: int
    ranking_resolution_bp: int
    sequence: Optional[str]
    genes: pd.DataFrame
    transcripts: pd.DataFrame
    exons: pd.DataFrame
    bins: pd.DataFrame


@dataclass
class NormalizedPrediction:
    interval_tracks: dict[str, np.ndarray]
    contact_map: Optional[np.ndarray]
    track_resolution_bp: dict[str, int]
    provenance: dict[str, Any]


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def default_fasta_path(organism: str) -> Optional[str]:
    organism = organism.lower()
    root = _project_root()
    if organism == "mouse":
        path = root / "data" / "raw" / "mm10.fa"
    else:
        path = root / "data" / "references" / "hg38.fa"
    return str(path) if path.exists() else None


def default_local_cooler_filename(organism: str) -> Optional[str]:
    if organism.lower() == "mouse":
        root = _project_root() / "data" / "raw" / "mouse_microc.mcool"
        return root.name if root.exists() else None
    return None


def _load_sequence_from_fasta(fasta_path: str, chrom: str, start: int, end: int) -> str:
    from pyfaidx import Fasta

    fa = Fasta(fasta_path)
    candidates = [chrom]
    if chrom.startswith("chr"):
        candidates.append(chrom[3:])
    else:
        candidates.append(f"chr{chrom}")
    for key in candidates:
        if key in fa:
            return str(fa[key][start:end]).upper()
    raise KeyError(f"Chromosome {chrom!r} not found in {fasta_path}")


def _fetch_ensembl_overlap(organism: str, chrom: str, start: int, end: int) -> pd.DataFrame:
    species = ENSEMBL_SPECIES.get(organism.lower(), organism.lower())
    region = f"{chrom}:{start + 1}-{end}"
    params = urllib.parse.urlencode([("feature", "gene"), ("feature", "transcript"), ("feature", "exon")])
    url = f"https://rest.ensembl.org/overlap/region/{species}/{region}?{params}"
    req = urllib.request.Request(url, headers={"Content-Type": "application/json", "Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=30) as resp:
        payload = json.loads(resp.read().decode())
    if not payload:
        return pd.DataFrame(columns=["feature_type", "id", "Parent", "start", "end", "strand", "external_name"])
    rows = []
    for rec in payload:
        rows.append(
            {
                "feature_type": rec.get("feature_type"),
                "id": rec.get("id"),
                "Parent": rec.get("Parent"),
                "start": int(rec.get("start", 0)) - 1,
                "end": int(rec.get("end", 0)),
                "strand": "+" if int(rec.get("strand", 1)) >= 0 else "-",
                "external_name": rec.get("external_name") or rec.get("gene_id") or rec.get("id"),
            }
        )
    return pd.DataFrame(rows)


def _normalize_annotations(annotation_df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    if annotation_df.empty:
        empty_gene = pd.DataFrame(columns=["gene_id", "gene_name", "chrom", "start", "end", "strand"])
        empty_tx = pd.DataFrame(columns=["transcript_id", "gene_id", "gene_name", "chrom", "start", "end", "strand", "tss"])
        empty_exon = pd.DataFrame(columns=["transcript_id", "gene_id", "chrom", "start", "end", "strand"])
        return empty_gene, empty_tx, empty_exon

    df = annotation_df.copy()
    if "Feature" in df.columns:
        feature_col = "Feature"
        gene_id_col = "gene_id"
        tx_id_col = "transcript_id"
        gene_name_col = "gene_name" if "gene_name" in df.columns else "gene_id"
        chrom_col = "seqname" if "seqname" in df.columns else "Chromosome"
        strand_col = "strand" if "strand" in df.columns else "Strand"
        start_col = "start" if "start" in df.columns else "Start"
        end_col = "end" if "end" in df.columns else "End"
    else:
        feature_col = "feature_type"
        gene_id_col = "id"
        tx_id_col = "id"
        gene_name_col = "external_name"
        chrom_col = "chrom"
        strand_col = "strand"
        start_col = "start"
        end_col = "end"

    genes_src = df[df[feature_col].isin(["gene", "Gene"])].copy()
    transcripts_src = df[df[feature_col].isin(["transcript", "Transcript"])].copy()
    exons_src = df[df[feature_col].isin(["exon", "Exon"])].copy()

    genes = pd.DataFrame(
        {
            "gene_id": genes_src.get(gene_id_col, pd.Series(dtype=str)).fillna("unknown"),
            "gene_name": genes_src.get(gene_name_col, pd.Series(dtype=str)).fillna("unknown"),
            "chrom": genes_src.get(chrom_col, pd.Series(dtype=str)),
            "start": genes_src.get(start_col, pd.Series(dtype=int)).astype(int),
            "end": genes_src.get(end_col, pd.Series(dtype=int)).astype(int),
            "strand": genes_src.get(strand_col, pd.Series(dtype=str)).fillna("+"),
        }
    )
    transcripts = pd.DataFrame(
        {
            "transcript_id": transcripts_src.get(tx_id_col, pd.Series(dtype=str)).fillna("unknown_tx"),
            "gene_id": transcripts_src.get("Parent", transcripts_src.get(gene_id_col, pd.Series(dtype=str))).fillna("unknown"),
            "gene_name": transcripts_src.get(gene_name_col, pd.Series(dtype=str)).fillna("unknown"),
            "chrom": transcripts_src.get(chrom_col, pd.Series(dtype=str)),
            "start": transcripts_src.get(start_col, pd.Series(dtype=int)).astype(int),
            "end": transcripts_src.get(end_col, pd.Series(dtype=int)).astype(int),
            "strand": transcripts_src.get(strand_col, pd.Series(dtype=str)).fillna("+"),
        }
    )
    transcripts["tss"] = np.where(transcripts["strand"] == "-", transcripts["end"], transcripts["start"])
    exons = pd.DataFrame(
        {
            "transcript_id": exons_src.get("Parent", exons_src.get(tx_id_col, pd.Series(dtype=str))).fillna("unknown_tx"),
            "gene_id": exons_src.get(gene_id_col, exons_src.get("Parent", pd.Series(dtype=str))).fillna("unknown"),
            "chrom": exons_src.get(chrom_col, pd.Series(dtype=str)),
            "start": exons_src.get(start_col, pd.Series(dtype=int)).astype(int),
            "end": exons_src.get(end_col, pd.Series(dtype=int)).astype(int),
            "strand": exons_src.get(strand_col, pd.Series(dtype=str)).fillna("+"),
        }
    )
    if not transcripts.empty and not genes.empty:
        gene_lookup = genes.drop_duplicates("gene_id").set_index("gene_id")
        transcripts["gene_name"] = transcripts["gene_id"].map(gene_lookup["gene_name"]).fillna(transcripts["gene_name"])
    return genes.drop_duplicates(), transcripts.drop_duplicates(), exons.drop_duplicates()


def _build_bins(chrom: str, start: int, end: int, resolution_bp: int) -> pd.DataFrame:
    edges = np.arange(start, end, resolution_bp, dtype=int)
    if len(edges) == 0 or edges[-1] != end:
        edges = np.append(edges, end)
    bins = pd.DataFrame({"chrom": chrom, "start": edges[:-1], "end": edges[1:]})
    bins["bin_index"] = np.arange(len(bins))
    return bins


def build_locus_context(
    request: LocusRequest,
    *,
    annotation_df: Optional[pd.DataFrame] = None,
) -> LocusContext:
    if request.sequence and not (request.anchor_interval or request.interval):
        raise ValueError("Synthetic sequence input requires an anchor interval for interpretation.")
    interval_str = request.anchor_interval or request.interval
    if not interval_str:
        raise ValueError("An interval or anchor interval is required.")

    chrom, input_start, input_end = parse_coordinates(interval_str)
    center = (input_start + input_end) // 2
    window_start = max(0, center - request.window_bp // 2)
    window_end = window_start + request.window_bp
    query_start, query_end = (parse_coordinates(request.interval)[1:3] if request.interval else (input_start, input_end))

    fasta_path = request.fasta_path or default_fasta_path(request.organism)
    seq = request.sequence
    if seq is None and fasta_path:
        seq = _load_sequence_from_fasta(fasta_path, chrom, window_start, window_end)

    if annotation_df is None:
        try:
            annotation_df = _fetch_ensembl_overlap(request.organism, chrom, window_start, window_end)
        except (urllib.error.URLError, TimeoutError, json.JSONDecodeError):
            annotation_df = pd.DataFrame()

    genes, transcripts, exons = _normalize_annotations(annotation_df)
    bins = _build_bins(chrom, window_start, window_end, request.ranking_resolution_bp)
    return LocusContext(
        request=request,
        chrom=chrom,
        interval_start=input_start,
        interval_end=input_end,
        window_start=window_start,
        window_end=window_end,
        query_start=query_start,
        query_end=query_end,
        ranking_resolution_bp=request.ranking_resolution_bp,
        sequence=seq,
        genes=genes,
        transcripts=transcripts,
        exons=exons,
        bins=bins,
    )


def _safe_array(values: Any) -> Optional[np.ndarray]:
    if values is None:
        return None
    if isinstance(values, np.ndarray):
        return values
    if hasattr(values, "values"):
        return np.asarray(values.values)
    try:
        return np.asarray(values)
    except Exception:
        return None


def _mean_interval_track(track_obj: Any) -> tuple[Optional[np.ndarray], Optional[int]]:
    arr = _safe_array(track_obj)
    if arr is None:
        return None, None
    if arr.ndim == 1:
        track = arr.astype(float)
    elif arr.ndim == 2:
        track = arr.mean(axis=1).astype(float)
    elif arr.ndim == 3:
        track = arr.mean(axis=(1, 2)).astype(float)
    else:
        return None, None
    resolution = getattr(track_obj, "resolution", None)
    return track, int(resolution) if resolution is not None else None


def _extract_ctcf_track(chip_tf_obj: Any) -> tuple[Optional[np.ndarray], Optional[int]]:
    arr = _safe_array(chip_tf_obj)
    if arr is None:
        return None, None
    if arr.ndim == 1:
        return arr.astype(float), getattr(chip_tf_obj, "resolution", None)

    metadata = getattr(chip_tf_obj, "metadata", None)
    if metadata is not None and hasattr(metadata, "__getitem__") and "transcription_factor" in metadata:
        mask = np.asarray(metadata["transcription_factor"] == "CTCF")
        if mask.any():
            if arr.ndim == 2:
                track = arr[:, mask].mean(axis=1).astype(float)
            else:
                track = arr.mean(axis=0)[:, mask].mean(axis=1).astype(float)
            return track, getattr(chip_tf_obj, "resolution", None)
    if arr.ndim == 2:
        return arr.mean(axis=1).astype(float), getattr(chip_tf_obj, "resolution", None)
    return arr.mean(axis=(1, 2)).astype(float), getattr(chip_tf_obj, "resolution", None)


def _extract_contact_map(contact_obj: Any) -> tuple[Optional[np.ndarray], Optional[int]]:
    arr = _safe_array(contact_obj)
    if arr is None:
        return None, None
    if arr.ndim == 2:
        mat = arr.astype(float)
    elif arr.ndim == 3:
        mat = arr.mean(axis=0).astype(float)
    else:
        return None, None
    return mat, getattr(contact_obj, "resolution", None)


class AlphaGenomePredictionAdapter:
    REQUESTED_OUTPUTS = [
        "RNA_SEQ",
        "CAGE",
        "ATAC",
        "DNASE",
        "CHIP_TF",
        "CONTACT_MAPS",
        "SPLICE_SITES",
        "SPLICE_SITE_USAGE",
    ]

    def __init__(self, connector: Optional[Any] = None):
        self.connector = connector or AlphaGenomeConnector()

    def predict_reference(self, context: LocusContext) -> NormalizedPrediction:
        if context.request.sequence is not None and context.sequence is not None and hasattr(self.connector, "predict_sequence"):
            raw = self.connector.predict_sequence(
                sequence=context.sequence,
                organism=context.request.organism.upper(),
                requested_outputs=self.REQUESTED_OUTPUTS,
                cell_type=context.request.cell_type,
            )
        else:
            raw = self.connector.predict_interval(
                chrom=context.chrom,
                start=context.window_start,
                end=context.window_end,
                organism=context.request.organism.upper(),
                requested_outputs=self.REQUESTED_OUTPUTS,
                cell_type=context.request.cell_type,
            )
        return self._normalize(raw, context.request)

    def predict_sequence(self, sequence: str, request: LocusRequest) -> NormalizedPrediction:
        raw = self.connector.predict_sequence(
            sequence=sequence,
            organism=request.organism.upper(),
            requested_outputs=self.REQUESTED_OUTPUTS,
            cell_type=request.cell_type,
        )
        return self._normalize(raw, request)

    def _normalize(self, raw: Any, request: LocusRequest) -> NormalizedPrediction:
        interval_tracks: dict[str, np.ndarray] = {}
        track_resolution_bp: dict[str, int] = {}
        provenance = {
            "backend": "AlphaGenome",
            "organism": request.organism,
            "cell_type": request.cell_type,
            "requested_outputs": list(self.REQUESTED_OUTPUTS),
            "available_modalities": [],
            "missing_modalities": [],
            "ontology_fallback_used": bool(getattr(self.connector, "last_ontology_fallback_used", False)),
        }

        for key, attr in [("rna_seq", "rna_seq"), ("cage", "cage"), ("atac", "atac"), ("dnase", "dnase"),
                          ("splice_sites", "splice_sites"), ("splice_site_usage", "splice_site_usage")]:
            obj = getattr(raw, attr, None)
            track, resolution = _mean_interval_track(obj)
            if track is not None:
                interval_tracks[key] = track
                track_resolution_bp[key] = int(resolution or request.ranking_resolution_bp)
                provenance["available_modalities"].append(key)
            else:
                provenance["missing_modalities"].append(key)

        ctcf_track, ctcf_resolution = _extract_ctcf_track(getattr(raw, "chip_tf", None))
        if ctcf_track is not None:
            interval_tracks["ctcf"] = ctcf_track
            track_resolution_bp["ctcf"] = int(ctcf_resolution or request.ranking_resolution_bp)
            provenance["available_modalities"].append("ctcf")
        else:
            provenance["missing_modalities"].append("ctcf")

        contact_map, contact_resolution = _extract_contact_map(getattr(raw, "contact_maps", None))
        if contact_map is not None:
            provenance["available_modalities"].append("contact_maps")
            track_resolution_bp["contact_maps"] = int(contact_resolution or request.ranking_resolution_bp)
        else:
            provenance["missing_modalities"].append("contact_maps")

        return NormalizedPrediction(
            interval_tracks=interval_tracks,
            contact_map=contact_map,
            track_resolution_bp=track_resolution_bp,
            provenance=provenance,
        )


def _resample_track(values: Optional[np.ndarray], source_resolution_bp: Optional[int], bins: pd.DataFrame) -> np.ndarray:
    if values is None or len(values) == 0:
        return np.zeros(len(bins), dtype=float)
    if source_resolution_bp is None or source_resolution_bp <= 0:
        source_resolution_bp = (bins.iloc[0]["end"] - bins.iloc[0]["start"])
    x_src = np.arange(len(values)) * source_resolution_bp
    target_resolution = int(bins.iloc[0]["end"] - bins.iloc[0]["start"])
    x_tgt = np.arange(len(bins)) * target_resolution
    return np.interp(x_tgt, x_src, values, left=float(values[0]), right=float(values[-1]))


def _normalize_vector(values: np.ndarray) -> np.ndarray:
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return arr
    finite = arr[np.isfinite(arr)]
    if finite.size == 0:
        return np.zeros_like(arr)
    lo = np.nanmin(finite)
    hi = np.nanmax(finite)
    if math.isclose(lo, hi):
        return np.zeros_like(arr)
    return np.clip((arr - lo) / (hi - lo), 0.0, 1.0)


def _contact_score_by_bin(contact_map: Optional[np.ndarray]) -> np.ndarray:
    if contact_map is None:
        return np.array([])
    n = contact_map.shape[0]
    scores = np.zeros(n, dtype=float)
    for i in range(n):
        row = np.asarray(contact_map[i], dtype=float).copy()
        lo = max(0, i - 1)
        hi = min(n, i + 2)
        row[lo:hi] = 0.0
        scores[i] = row.sum()
    return _normalize_vector(scores)


def _insulation_track(contact_map: Optional[np.ndarray]) -> np.ndarray:
    if contact_map is None:
        return np.array([])
    baseline = np.zeros_like(contact_map)
    delta = compute_insulation_delta(baseline, contact_map, window_bins=max(2, min(10, contact_map.shape[0] // 10)))
    return np.nan_to_num(np.abs(delta), nan=0.0)


def _domain_ids(boundary_track: np.ndarray, threshold_quantile: float = 0.85) -> np.ndarray:
    if boundary_track.size == 0:
        return np.zeros(0, dtype=int)
    thresh = float(np.quantile(boundary_track, threshold_quantile))
    domain = np.zeros(len(boundary_track), dtype=int)
    current = 0
    for i, score in enumerate(boundary_track):
        if i > 0 and score >= thresh:
            current += 1
        domain[i] = current
    return domain


def _match_cell_type(cell_type: str, organism: str) -> tuple[Optional[dict[str, Any]], Optional[dict[str, Any]]]:
    root = _project_root() / "notebooks"
    prefix = "mouse" if organism.lower() == "mouse" else "human"
    tracks_path = root / f"{prefix}_tracks_metadata.csv"
    contacts_path = root / f"{prefix}_contact_maps_metadata.csv"
    result_track = None
    result_contact = None
    query = (cell_type or "").strip().lower()

    def _search(path: Path) -> Optional[dict[str, Any]]:
        if not path.exists():
            return None
        df = pd.read_csv(path)
        if df.empty:
            return None
        text = (
            df.get("biosample_name", pd.Series(dtype=str)).fillna("").astype(str).str.lower()
            + " "
            + df.get("ontology_curie", pd.Series(dtype=str)).fillna("").astype(str).str.lower()
            + " "
            + df.get("name", pd.Series(dtype=str)).fillna("").astype(str).str.lower()
        )
        exact = df[text.str.contains(query, regex=False)]
        if not exact.empty:
            row = exact.iloc[0].to_dict()
            row["match_type"] = "exact"
            return row
        if not df.empty:
            row = df.iloc[0].to_dict()
            row["match_type"] = "fallback"
            return row
        return None

    result_track = _search(tracks_path)
    result_contact = _search(contacts_path)
    return result_track, result_contact


def _observed_contact_support(context: LocusContext) -> tuple[np.ndarray, dict[str, Any]]:
    filename = context.request.local_cooler_filename or default_local_cooler_filename(context.request.organism)
    if not filename:
        return np.zeros(len(context.bins), dtype=float), {"source": None, "status": "missing"}
    try:
        clr = load_cooler(filename, resolution=context.request.local_cooler_resolution)
        matrix = clr.matrix(balance=True).fetch(f"{context.chrom}:{context.window_start}-{context.window_end}")
        scores = _contact_score_by_bin(np.nan_to_num(matrix, nan=0.0))
        if scores.size != len(context.bins):
            scores = np.interp(np.arange(len(context.bins)), np.linspace(0, len(context.bins) - 1, len(scores)), scores)
        return scores, {
            "source": filename,
            "status": "loaded",
            "resolution_bp": context.request.local_cooler_resolution,
        }
    except Exception as exc:
        return np.zeros(len(context.bins), dtype=float), {
            "source": filename,
            "status": "error",
            "error": str(exc),
        }


def load_external_evidence(context: LocusContext) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    track_match, contact_match = _match_cell_type(context.request.cell_type, context.request.organism)
    observed_contact, contact_info = _observed_contact_support(context)
    evidence = {"observed_contact": observed_contact}
    provenance = {
        "track_metadata_match": track_match,
        "contact_metadata_match": contact_match,
        "observed_contact": contact_info,
        "biosample_fallbacks": [],
    }
    for match in (track_match, contact_match):
        if match and match.get("match_type") == "fallback":
            provenance["biosample_fallbacks"].append(match.get("biosample_name") or match.get("name"))
    return evidence, provenance


def _promoter_bins(context: LocusContext) -> pd.DataFrame:
    if context.transcripts.empty:
        return pd.DataFrame(columns=["transcript_id", "gene_id", "gene_name", "tss", "bin_index", "domain_id"])
    promoters = context.transcripts.copy()
    promoters["bin_index"] = ((promoters["tss"] - context.window_start) // context.ranking_resolution_bp).clip(lower=0)
    promoters["bin_index"] = promoters["bin_index"].astype(int)
    promoters = promoters[(promoters["bin_index"] >= 0) & (promoters["bin_index"] < len(context.bins))]
    return promoters


def _distance_score(a: int, b: int) -> float:
    return 1.0 / (1.0 + abs(a - b))


def _build_gene_summaries(
    context: LocusContext,
    bins_df: pd.DataFrame,
    promoters: pd.DataFrame,
    prediction: NormalizedPrediction,
    links_df: pd.DataFrame,
) -> list[GeneActivitySummary]:
    summaries: list[GeneActivitySummary] = []
    rna = bins_df["rna_seq"].values if "rna_seq" in bins_df else np.zeros(len(bins_df))
    splice = bins_df["splice_site_usage"].values if "splice_site_usage" in bins_df else np.zeros(len(bins_df))
    gene_lookup = context.genes.drop_duplicates("gene_id").set_index("gene_id") if not context.genes.empty else None
    for _, promoter in promoters.iterrows():
        gene_id = promoter["gene_id"]
        gene_name = promoter.get("gene_name") or gene_id
        gene_row = gene_lookup.loc[gene_id] if gene_lookup is not None and gene_id in gene_lookup.index else None
        promoter_bin = int(promoter["bin_index"])
        local_links = links_df[links_df["gene_id"] == gene_id].sort_values("confidence", ascending=False)
        nearest = np.argsort(np.abs(bins_df["bin_index"].values - promoter_bin))[:3].tolist()
        summaries.append(
            GeneActivitySummary(
                gene_id=str(gene_id),
                gene_symbol=str(gene_name),
                chrom=context.chrom if gene_row is None else str(gene_row["chrom"]),
                start=int(context.window_start if gene_row is None else gene_row["start"]),
                end=int(context.window_end if gene_row is None else gene_row["end"]),
                domain_id=int(bins_df.loc[promoter_bin, "3d_domain_id"]) if promoter_bin < len(bins_df) else 0,
                promoter_bin=promoter_bin,
                promoter_activity_score=float(bins_df.loc[promoter_bin, "promoter_activity_score"]),
                nearest_elements=nearest,
                supporting_elements=local_links["enhancer_bin"].astype(int).tolist(),
                linked_element_bins=local_links["enhancer_bin"].astype(int).tolist(),
                predicted_expression_change=float(rna[promoter_bin]) if promoter_bin < len(rna) else 0.0,
                predicted_splice_change=float(splice[promoter_bin]) if promoter_bin < len(splice) else 0.0,
                linked_element_score=float(local_links["confidence"].sum()),
                evidence_sources=sorted(set(local_links["evidence_source"].tolist() + ["alphagenome"])),
            )
        )
    return summaries


def _classify_element(promoter_score: float, sequence_score: float, ctcf_score: float) -> str:
    if promoter_score >= 0.6:
        return "promoter"
    if ctcf_score >= 0.7 and sequence_score < 0.35:
        return "structural"
    if sequence_score >= 0.45:
        return "enhancer_candidate"
    return "inert"


def _serialize_df(df: pd.DataFrame) -> list[dict[str, Any]]:
    if df.empty:
        return []
    return json.loads(df.to_json(orient="records"))


def rank_regulatory_activity(
    context: LocusContext,
    prediction: NormalizedPrediction,
    external_evidence: Optional[dict[str, np.ndarray]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[GeneActivitySummary], dict[str, Any]]:
    external_evidence = external_evidence or {}
    bins = context.bins.copy()
    for name in ["rna_seq", "cage", "atac", "dnase", "splice_sites", "splice_site_usage", "ctcf"]:
        track = prediction.interval_tracks.get(name)
        bins[name] = _resample_track(track, prediction.track_resolution_bp.get(name), bins)
        bins[name] = _normalize_vector(bins[name].values)

    bins["observed_contact_score"] = _normalize_vector(external_evidence.get("observed_contact", np.zeros(len(bins))))
    bins["contact_score"] = _normalize_vector(_resample_track(_contact_score_by_bin(prediction.contact_map), None, bins))
    bins["boundary_context_score"] = _normalize_vector(_resample_track(_insulation_track(prediction.contact_map), None, bins))
    bins["ctcf_score"] = bins["ctcf"]
    bins["sequence_activity_score"] = _normalize_vector(
        0.35 * bins["atac"].values
        + 0.20 * bins["dnase"].values
        + 0.25 * bins["cage"].values
        + 0.20 * bins["rna_seq"].values
    )
    bins["promoter_activity_score"] = _normalize_vector(0.7 * bins["cage"].values + 0.3 * bins["rna_seq"].values)
    bins["cell_type_support_score"] = _normalize_vector(0.6 * bins["observed_contact_score"].values + 0.4 * bins["atac"].values)
    bins["3d_domain_id"] = _domain_ids(bins["boundary_context_score"].values)

    promoters = _promoter_bins(context)
    promoter_bin_set = set(promoters["bin_index"].astype(int).tolist())
    bins["is_promoter_bin"] = bins["bin_index"].isin(promoter_bin_set)

    links = []
    for _, row in bins.iterrows():
        if row["bin_index"] in promoter_bin_set:
            continue
        candidate_type = _classify_element(float(row["promoter_activity_score"]), float(row["sequence_activity_score"]), float(row["ctcf_score"]))
        if candidate_type == "inert":
            continue
        domain_promoters = promoters[promoters["bin_index"].map(lambda idx: bins.loc[int(idx), "3d_domain_id"] if int(idx) < len(bins) else -1) == row["3d_domain_id"]]
        if domain_promoters.empty:
            domain_promoters = promoters
        for _, promoter in domain_promoters.iterrows():
            promoter_bin = int(promoter["bin_index"])
            contact = float(row["contact_score"])
            confidence = float(row["sequence_activity_score"] * max(contact, _distance_score(int(row["bin_index"]), promoter_bin)))
            links.append(
                {
                    "enhancer_bin": int(row["bin_index"]),
                    "promoter_bin": promoter_bin,
                    "gene_id": promoter["gene_id"],
                    "gene_name": promoter.get("gene_name") or promoter["gene_id"],
                    "confidence": confidence,
                    "contact_score": contact,
                    "evidence_source": "activity_x_contact",
                }
            )
    links_df = pd.DataFrame(links)
    if not links_df.empty:
        top_links = links_df.sort_values("confidence", ascending=False).groupby("enhancer_bin").head(1)
        promoter_link_lookup = top_links.set_index("enhancer_bin")
    else:
        top_links = pd.DataFrame(columns=["enhancer_bin", "gene_id", "gene_name", "confidence", "evidence_source"])
        promoter_link_lookup = top_links.set_index("enhancer_bin")

    bins["promoter_link_score"] = bins["bin_index"].map(lambda idx: float(promoter_link_lookup.loc[idx, "confidence"]) if idx in promoter_link_lookup.index else 0.0)
    bins["element_type"] = bins.apply(
        lambda row: "promoter" if row["bin_index"] in promoter_bin_set else _classify_element(float(row["promoter_activity_score"]), float(row["sequence_activity_score"]), float(row["ctcf_score"])),
        axis=1,
    )
    bins["activity_score"] = (
        0.35 * bins["sequence_activity_score"]
        + 0.25 * bins["promoter_link_score"]
        + 0.20 * bins["contact_score"]
        + 0.10 * bins["cell_type_support_score"]
        + 0.10 * bins["boundary_context_score"]
    )
    structural_mask = (bins["element_type"] == "structural") & (bins["sequence_activity_score"] < 0.35)
    bins.loc[structural_mask, "activity_score"] *= 0.5

    linked_gene_ids = []
    linked_gene_symbols = []
    evidence_sources = []
    for idx in bins["bin_index"]:
        if idx in promoter_link_lookup.index:
            linked_gene_ids.append([str(promoter_link_lookup.loc[idx, "gene_id"])])
            linked_gene_symbols.append([str(promoter_link_lookup.loc[idx, "gene_name"])])
            evidence_sources.append([str(promoter_link_lookup.loc[idx, "evidence_source"]), "alphagenome"])
        else:
            linked_gene_ids.append([])
            linked_gene_symbols.append([])
            evidence_sources.append(["alphagenome"])
    bins["linked_gene_ids"] = linked_gene_ids
    bins["linked_gene_symbols"] = linked_gene_symbols
    bins["evidence_sources"] = evidence_sources

    gene_summaries = _build_gene_summaries(context, bins, promoters, prediction, links_df)
    summary = {
        "n_bins": int(len(bins)),
        "n_promoters": int(len(promoters)),
        "n_links": int(len(links_df)),
        "n_candidate_elements": int((bins["element_type"] != "inert").sum()),
    }
    return bins.sort_values("activity_score", ascending=False).reset_index(drop=True), top_links.sort_values("confidence", ascending=False).reset_index(drop=True), gene_summaries, summary


def _apply_edit_to_sequence(context: LocusContext, sequence: str, edit: dict[str, Any]) -> str:
    if context.sequence is None:
        raise ValueError("Perturbation mode requires a reference sequence.")
    mode = str(edit.get("mode", "")).lower()
    start = int(edit["start"])
    end = int(edit.get("end", start))
    insert_seq = str(edit.get("sequence", "")).upper()
    rel_start = start - context.window_start
    rel_end = end - context.window_start
    if rel_start < 0 or rel_end > len(sequence):
        raise ValueError("Edit lies outside the modeled sequence window.")

    if mode == "deletion":
        edited = sequence[:rel_start] + sequence[rel_end:]
        fill = sequence[-(rel_end - rel_start):] if rel_end > rel_start else ""
        edited = edited + fill
    elif mode == "substitution":
        if end <= start:
            raise ValueError("Substitution requires end > start.")
        edited = sequence[:rel_start] + insert_seq + sequence[rel_end:]
    elif mode == "insertion":
        edited = sequence[:rel_start] + insert_seq + sequence[rel_start:]
        edited = edited[:len(sequence)]
    else:
        raise ValueError(f"Unsupported edit mode: {mode}")

    if len(edited) < len(sequence):
        edited = edited + sequence[len(edited):len(sequence)]
    return edited[:len(sequence)]


def compare_reference_vs_edit(
    context: LocusContext,
    adapter: AlphaGenomePredictionAdapter,
    reference_bins: pd.DataFrame,
    edit: dict[str, Any],
    external_evidence: Optional[dict[str, np.ndarray]] = None,
) -> tuple[pd.DataFrame, pd.DataFrame, list[GeneActivitySummary], dict[str, Any], dict[str, Any], NormalizedPrediction]:
    edited_sequence = _apply_edit_to_sequence(context, context.sequence or "", edit)
    edited_prediction = adapter.predict_sequence(edited_sequence, context.request)
    edited_bins, edited_links, edited_gene_summaries, edited_summary = rank_regulatory_activity(context, edited_prediction, external_evidence)
    merged = reference_bins[["bin_index", "activity_score"]].merge(
        edited_bins[["bin_index", "activity_score"]],
        on="bin_index",
        suffixes=("_ref", "_edit"),
    )
    delta_lookup = merged.set_index("bin_index")
    edited_bins["delta_activity_score"] = edited_bins["bin_index"].map(
        lambda idx: float(delta_lookup.loc[idx, "activity_score_edit"] - delta_lookup.loc[idx, "activity_score_ref"]) if idx in delta_lookup.index else 0.0
    )
    diff_summary = {
        "mode": edit.get("mode"),
        "gained_bins": edited_bins.sort_values("delta_activity_score", ascending=False).head(5)["bin_index"].astype(int).tolist(),
        "lost_bins": edited_bins.sort_values("delta_activity_score", ascending=True).head(5)["bin_index"].astype(int).tolist(),
    }
    return edited_bins, edited_links, edited_gene_summaries, edited_summary, diff_summary, edited_prediction


def _candidate_from_row(row: pd.Series) -> CandidateElement:
    return CandidateElement(
        chrom=str(row["chrom"]),
        start=int(row["start"]),
        end=int(row["end"]),
        bin_index=int(row["bin_index"]),
        activity_score=float(row["activity_score"]),
        sequence_activity_score=float(row["sequence_activity_score"]),
        promoter_link_score=float(row["promoter_link_score"]),
        contact_score=float(row["contact_score"]),
        boundary_context_score=float(row["boundary_context_score"]),
        cell_type_support_score=float(row["cell_type_support_score"]),
        promoter_activity_score=float(row["promoter_activity_score"]),
        ctcf_score=float(row["ctcf_score"]),
        observed_contact_score=float(row["observed_contact_score"]),
        d3_domain_id=int(row["3d_domain_id"]),
        element_type=str(row["element_type"]),
        linked_gene_ids=list(row["linked_gene_ids"]),
        linked_gene_symbols=list(row["linked_gene_symbols"]),
        evidence_sources=list(row["evidence_sources"]),
        delta_activity_score=(None if pd.isna(row.get("delta_activity_score")) else float(row["delta_activity_score"])),
    )


def build_manifest(
    context: LocusContext,
    reference_bins: pd.DataFrame,
    promoter_links: pd.DataFrame,
    gene_summaries: list[GeneActivitySummary],
    prediction: NormalizedPrediction,
    external_provenance: dict[str, Any],
    scoring_summary: dict[str, Any],
    *,
    edited_bins: Optional[pd.DataFrame] = None,
    edited_links: Optional[pd.DataFrame] = None,
    edited_gene_summaries: Optional[list[GeneActivitySummary]] = None,
    edited_prediction: Optional[NormalizedPrediction] = None,
    perturbation_summary: Optional[dict[str, Any]] = None,
) -> dict[str, Any]:
    manifest = {
        "input": {
            "organism": context.request.organism,
            "assembly": context.request.assembly,
            "cell_type": context.request.cell_type,
            "interval": context.request.interval,
            "anchor_interval": context.request.anchor_interval,
            "source_mode": "synthetic_sequence" if context.request.sequence else "reference_interval",
            "window": f"{context.chrom}:{context.window_start}-{context.window_end}",
        },
        "candidate_elements": [asdict(_candidate_from_row(row)) for _, row in reference_bins.iterrows()],
        "promoter_links": _serialize_df(promoter_links),
        "gene_summaries": [asdict(item) for item in gene_summaries],
        "model_provenance": {
            "predicted": prediction.provenance,
            "external": external_provenance,
        },
        "summary": scoring_summary,
    }
    if edited_bins is not None:
        manifest["perturbation"] = {
            "summary": perturbation_summary or {},
            "candidate_elements": [asdict(_candidate_from_row(row)) for _, row in edited_bins.iterrows()],
            "promoter_links": _serialize_df(edited_links if edited_links is not None else pd.DataFrame()),
            "gene_summaries": [asdict(item) for item in (edited_gene_summaries or [])],
            "model_provenance": edited_prediction.provenance if edited_prediction is not None else {},
        }
    return manifest


def _slugify(text: str) -> str:
    keep = [c.lower() if c.isalnum() else "_" for c in text]
    slug = "".join(keep)
    while "__" in slug:
        slug = slug.replace("__", "_")
    return slug.strip("_")


def manifest_slug(context: LocusContext) -> str:
    label = f"{context.request.organism}_{context.request.cell_type}_{context.chrom}_{context.query_start}_{context.query_end}"
    return _slugify(label)


def save_outputs(
    context: LocusContext,
    manifest: dict[str, Any],
    reference_bins: pd.DataFrame,
    gene_summaries: list[GeneActivitySummary],
    *,
    edited_bins: Optional[pd.DataFrame] = None,
) -> dict[str, str]:
    slug = manifest_slug(context)
    processed_dir = _project_root() / "data" / "processed"
    media_dir = _project_root() / "media"
    processed_dir.mkdir(parents=True, exist_ok=True)
    media_dir.mkdir(parents=True, exist_ok=True)

    manifest_path = processed_dir / f"regulatory_activity_{slug}.json"
    with open(manifest_path, "w") as f:
        json.dump(manifest, f, indent=2)

    fig_paths = {}
    overview_path = media_dir / f"regulatory_activity_{slug}_overview.png"
    _plot_overview(context, reference_bins, gene_summaries, overview_path)
    fig_paths["overview"] = str(overview_path)

    ranking_path = media_dir / f"regulatory_activity_{slug}_top_candidates.png"
    _plot_top_candidates(reference_bins, ranking_path)
    fig_paths["top_candidates"] = str(ranking_path)

    if edited_bins is not None:
        delta_path = media_dir / f"regulatory_activity_{slug}_delta.png"
        _plot_delta(edited_bins, delta_path)
        fig_paths["delta"] = str(delta_path)

    fig_paths["manifest"] = str(manifest_path)
    return fig_paths


def _plot_overview(context: LocusContext, reference_bins: pd.DataFrame, gene_summaries: list[GeneActivitySummary], out_path: Path) -> None:
    fig, axes = plt.subplots(3, 1, figsize=(14, 8), sharex=True)
    x = reference_bins["start"].values
    axes[0].fill_between(x, reference_bins["activity_score"].values, color="#1f77b4", alpha=0.6)
    axes[0].plot(x, reference_bins["sequence_activity_score"].values, color="#d62728", linewidth=1.0, label="sequence")
    axes[0].plot(x, reference_bins["contact_score"].values, color="#2ca02c", linewidth=1.0, label="contact")
    axes[0].set_ylabel("Score")
    axes[0].set_title(
        textwrap.fill(
            f"Regulatory activity overview: {context.request.organism} {context.request.cell_type} {context.chrom}:{context.query_start}-{context.query_end}",
            100,
        )
    )
    axes[0].legend(loc="upper right", fontsize=8)

    axes[1].bar(x, reference_bins["promoter_link_score"].values, width=context.ranking_resolution_bp, color="#9467bd", alpha=0.7)
    axes[1].plot(x, reference_bins["boundary_context_score"].values, color="#8c564b", linewidth=1.0, label="boundary")
    axes[1].set_ylabel("Link / boundary")
    axes[1].legend(loc="upper right", fontsize=8)

    axes[2].scatter(x, reference_bins["activity_score"].values, c=reference_bins["3d_domain_id"].values, cmap="tab20", s=14)
    for gene in gene_summaries:
        axes[2].axvline(context.window_start + gene.promoter_bin * context.ranking_resolution_bp, color="black", alpha=0.15, linewidth=0.8)
        axes[2].text(
            context.window_start + gene.promoter_bin * context.ranking_resolution_bp,
            1.02,
            gene.gene_symbol,
            rotation=90,
            va="bottom",
            ha="center",
            fontsize=7,
        )
    axes[2].set_ylabel("Activity")
    axes[2].set_xlabel("Genomic position")
    axes[2].set_ylim(0, 1.1)

    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_top_candidates(reference_bins: pd.DataFrame, out_path: Path) -> None:
    top = reference_bins.head(15).copy().sort_values("activity_score")
    labels = [f"{row.element_type}:{int(row.bin_index)}" for row in top.itertuples()]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.barh(labels, top["activity_score"].values, color="#1f77b4", alpha=0.8)
    ax.set_xlabel("Activity score")
    ax.set_title("Top ranked regulatory sections")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def _plot_delta(edited_bins: pd.DataFrame, out_path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12, 4))
    ax.bar(edited_bins["start"].values, edited_bins["delta_activity_score"].values, width=(edited_bins["end"] - edited_bins["start"]).values, color="#ff7f0e", alpha=0.8)
    ax.axhline(0.0, color="black", linewidth=0.8)
    ax.set_xlabel("Genomic position")
    ax.set_ylabel("Delta activity")
    ax.set_title("Perturbation delta activity")
    fig.tight_layout()
    fig.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)

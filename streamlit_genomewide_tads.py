from __future__ import annotations

import math
from pathlib import Path
from typing import Optional
import os
import tomllib

import cooler
import numpy as np
import pandas as pd
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import streamlit as st


st.set_page_config(page_title="Genome-wide TAD Browser", layout="wide")

DEFAULT_PREFIX = Path("data/processed/genomewide_tads/mouse_microc_100kb")
DEFAULT_COOLER = Path("data/raw/mouse_microc.mcool")
DEFAULT_MANIFEST = Path("data/external/mm10/manifest.toml")


@st.cache_data
def load_callset(prefix_str: str) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, Path]:
    prefix = Path(prefix_str)
    tads = pd.read_csv(prefix.with_suffix(".tads.tsv"), sep="\t")
    boundaries = pd.read_csv(prefix.with_suffix(".boundaries.tsv"), sep="\t")
    tad_summary = pd.read_csv(prefix.with_suffix(".tad_summary.tsv"), sep="\t")
    return tads, boundaries, tad_summary, prefix


@st.cache_resource
def load_cooler(cooler_path: str, resolution: int) -> cooler.Cooler:
    return cooler.Cooler(f"{cooler_path}::resolutions/{resolution}")


@st.cache_data
def load_manifest(path_str: str = str(DEFAULT_MANIFEST)) -> dict:
    path = Path(path_str)
    if not path.exists():
        return {}
    with path.open("rb") as handle:
        return tomllib.load(handle)


def _optional_pybigwig():
    try:
        import pyBigWig  # type: ignore

        return pyBigWig
    except Exception:
        return None


def _safe_int(x: object) -> int:
    return int(float(x))


def _read_tabular_track(path: Path) -> pd.DataFrame:
    sep = "," if path.suffix.lower() == ".csv" else "\t"
    # Try headered first, then raw BED-style.
    try:
        df = pd.read_csv(path, sep=sep, comment="#")
    except Exception:
        df = pd.read_csv(path, sep=sep, comment="#", header=None)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "signal"])

    lower_cols = {str(c).lower(): c for c in df.columns}
    if {"chrom", "start", "end"}.issubset(lower_cols):
        chrom_col = lower_cols["chrom"]
        start_col = lower_cols["start"]
        end_col = lower_cols["end"]
        signal_candidates = [
            c for c in df.columns
            if str(c).lower() in {"signal", "score", "value", "tpm", "fpkm", "rpkm"}
        ]
        if not signal_candidates:
            signal_candidates = [c for c in df.columns if c not in {chrom_col, start_col, end_col, lower_cols.get("gene_name", "__none__")}]
        signal_col = signal_candidates[0] if signal_candidates else None
        out = pd.DataFrame(
            {
                "chrom": df[chrom_col].astype(str),
                "start": pd.to_numeric(df[start_col], errors="coerce"),
                "end": pd.to_numeric(df[end_col], errors="coerce"),
                "signal": pd.to_numeric(df[signal_col], errors="coerce") if signal_col is not None else 1.0,
            }
        )
        out = out.dropna(subset=["start", "end"])
        out["start"] = out["start"].astype(int)
        out["end"] = out["end"].astype(int)
        out["signal"] = pd.to_numeric(out["signal"], errors="coerce").fillna(0.0)
        return out

    df = pd.read_csv(path, sep=sep, comment="#", header=None)
    if df.shape[1] < 3:
        raise ValueError(f"Unsupported track format for {path}")

    chrom_col, start_col, end_col = 0, 1, 2
    signal_series = None
    if df.shape[1] == 3:
        signal_series = pd.Series(1.0, index=df.index)
    else:
        candidate_indices = list(range(3, min(df.shape[1], 10)))
        numeric_candidates = []
        for idx in candidate_indices:
            numeric = pd.to_numeric(df[idx], errors="coerce")
            numeric_candidates.append((idx, numeric.notna().sum(), numeric))
        numeric_candidates.sort(key=lambda x: x[1], reverse=True)
        if numeric_candidates and numeric_candidates[0][1] > 0:
            signal_series = numeric_candidates[0][2]
        else:
            signal_series = pd.Series(1.0, index=df.index)

    out = pd.DataFrame(
        {
            "chrom": df[chrom_col].astype(str),
            "start": pd.to_numeric(df[start_col], errors="coerce"),
            "end": pd.to_numeric(df[end_col], errors="coerce"),
            "signal": pd.to_numeric(signal_series, errors="coerce"),
        }
    )
    out = out.dropna(subset=["start", "end"])
    out["start"] = out["start"].astype(int)
    out["end"] = out["end"].astype(int)
    out["signal"] = out["signal"].fillna(0.0)
    return out


@st.cache_data
def load_signal_intervals(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in {".bw", ".bigwig", ".bigWig"}:
        pyBigWig = _optional_pybigwig()
        if pyBigWig is None:
            raise RuntimeError("pyBigWig is not installed in the current environment.")
        bw = pyBigWig.open(str(path))
        chroms = bw.chroms()
        rows: list[dict[str, object]] = []
        for chrom, chrom_len in chroms.items():
            intervals = bw.intervals(chrom)
            if not intervals:
                continue
            for start, end, value in intervals:
                rows.append({"chrom": chrom, "start": int(start), "end": int(end), "signal": float(value)})
        bw.close()
        return pd.DataFrame(rows, columns=["chrom", "start", "end", "signal"])

    return _read_tabular_track(path)


def _parse_gtf_attributes(attr_str: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for part in attr_str.strip().split(";"):
        part = part.strip()
        if not part:
            continue
        if " " not in part:
            continue
        key, value = part.split(" ", 1)
        result[key] = value.strip().strip('"')
    return result


@st.cache_data
def load_gene_annotations(path_str: str) -> pd.DataFrame:
    path = Path(path_str)
    if not path.exists():
        raise FileNotFoundError(path)

    if path.suffix.lower() in {".bed", ".bed12"}:
        df = pd.read_csv(path, sep="\t", comment="#", header=None)
        if df.shape[1] < 3:
            raise ValueError(f"Unsupported BED format in {path}")
        out = pd.DataFrame(
            {
                "chrom": df.iloc[:, 0].astype(str),
                "start": df.iloc[:, 1].astype(int),
                "end": df.iloc[:, 2].astype(int),
                "gene_name": df.iloc[:, 3].astype(str) if df.shape[1] > 3 else "feature",
                "strand": df.iloc[:, 5].astype(str) if df.shape[1] > 5 else ".",
            }
        )
        return out

    rows: list[dict[str, object]] = []
    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            attrs = _parse_gtf_attributes(parts[8])
            rows.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[3]) - 1,
                    "end": int(parts[4]),
                    "gene_name": attrs.get("gene_name", attrs.get("gene_id", "gene")),
                    "strand": parts[6],
                }
            )
    return pd.DataFrame(rows, columns=["chrom", "start", "end", "gene_name", "strand"])


def validate_track(path_str: str, expected_chr_prefix: str = "chr", kind: str = "signal") -> list[str]:
    path = Path(path_str)
    if not path.exists():
        return ["File does not exist"]

    issues: list[str] = []
    suffixes = {s.lower() for s in path.suffixes}

    if kind == "genes":
        try:
            genes = load_gene_annotations(path_str)
            if genes.empty:
                issues.append("No gene rows loaded")
            elif not genes["chrom"].astype(str).str.startswith(expected_chr_prefix).all():
                issues.append("Chromosome naming mismatch (expected 'chr*')")
        except Exception as exc:
            return [f"Failed to read: {exc}"]
        return issues

    if ".bigwig" in suffixes or ".bw" in suffixes:
        pyBigWig = _optional_pybigwig()
        if pyBigWig is None:
            issues.append("pyBigWig not installed")
            return issues
        try:
            bw = pyBigWig.open(str(path))
            chroms = list(bw.chroms().keys())
            bw.close()
            if not chroms:
                issues.append("No chromosomes found in bigWig")
            elif not all(str(c).startswith(expected_chr_prefix) for c in chroms[:20]):
                issues.append("Chromosome naming mismatch (expected 'chr*')")
        except Exception as exc:
            return [f"Failed to read: {exc}"]
        return issues

    try:
        df = pd.read_csv(path, sep="\t", comment="#", nrows=100)
        lower_cols = {str(c).lower() for c in df.columns}
        if {"chrom", "start", "end"}.issubset(lower_cols):
            chrom_col = next(c for c in df.columns if str(c).lower() == "chrom")
            if not df[chrom_col].astype(str).str.startswith(expected_chr_prefix).all():
                issues.append("Chromosome naming mismatch (expected 'chr*')")
            signal_cols = [c for c in df.columns if str(c).lower() not in {"chrom", "start", "end", "gene_name"}]
            if signal_cols and not pd.to_numeric(df[signal_cols[0]], errors="coerce").notnull().any():
                issues.append("No numeric signal/score column detected")
            return issues
    except Exception:
        pass

    try:
        df = pd.read_csv(path, sep="\t", comment="#", nrows=100, header=None)
    except Exception as exc:
        return [f"Failed to read: {exc}"]

    if df.empty:
        return ["File is empty"]
    if not df[0].astype(str).str.startswith(expected_chr_prefix).all():
        issues.append("Chromosome naming mismatch (expected 'chr*')")
    if df.shape[1] < 3:
        issues.append("Missing coordinate columns")
        return issues

    signal_ok = False
    for idx in range(3, min(df.shape[1], 10)):
        if pd.to_numeric(df[idx], errors="coerce").notnull().any():
            signal_ok = True
            break
    if not signal_ok:
        issues.append("No numeric signal/score column detected")
    return issues


def get_assay_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if any(col.startswith(prefix) for prefix in ("atac_", "chip_", "rna_"))]


def filter_intervals(df: pd.DataFrame, chrom: str, start: int, end: int) -> pd.DataFrame:
    if df.empty:
        return df.copy()
    return df[(df["chrom"] == chrom) & (df["end"] > start) & (df["start"] < end)].copy()


def interval_weighted_mean(df: pd.DataFrame, chrom: str, start: int, end: int) -> float:
    locus = filter_intervals(df, chrom, start, end)
    if locus.empty:
        return 0.0
    width = max(1, end - start)
    overlap_start = np.maximum(locus["start"].to_numpy(dtype=int), start)
    overlap_end = np.minimum(locus["end"].to_numpy(dtype=int), end)
    overlap = np.clip(overlap_end - overlap_start, 0, None)
    if not np.any(overlap > 0):
        return 0.0
    signal = locus["signal"].to_numpy(dtype=float)
    return float(np.sum(overlap * signal) / width)


def interval_max_signal(df: pd.DataFrame, chrom: str, start: int, end: int) -> float:
    locus = filter_intervals(df, chrom, start, end)
    if locus.empty:
        return 0.0
    return float(pd.to_numeric(locus["signal"], errors="coerce").fillna(0.0).max())


def interval_count(df: pd.DataFrame, chrom: str, start: int, end: int) -> int:
    return int(len(filter_intervals(df, chrom, start, end)))


def interval_coverage_fraction(df: pd.DataFrame, chrom: str, start: int, end: int) -> float:
    locus = filter_intervals(df, chrom, start, end)
    if locus.empty:
        return 0.0
    intervals = sorted((max(start, int(r.start)), min(end, int(r.end))) for r in locus.itertuples())
    merged: list[tuple[int, int]] = []
    for s, e in intervals:
        if e <= s:
            continue
        if not merged or s > merged[-1][1]:
            merged.append((s, e))
        else:
            merged[-1] = (merged[-1][0], max(merged[-1][1], e))
    covered = sum(e - s for s, e in merged)
    return float(covered / max(1, end - start))


def genes_with_rna_values(genes: pd.DataFrame, rna_df: Optional[pd.DataFrame], chrom: str) -> pd.DataFrame:
    if genes.empty:
        return genes.copy()
    out = genes.copy()
    values = []
    for gene in out.itertuples():
        if rna_df is None:
            values.append(0.0)
        else:
            values.append(interval_max_signal(rna_df, chrom, int(gene.start), int(gene.end)))
    out["rna_signal"] = values
    return out.sort_values("rna_signal", ascending=False)


def sample_signal_track(
    df: pd.DataFrame,
    chrom: str,
    start: int,
    end: int,
    n_bins: int = 500,
) -> pd.DataFrame:
    locus = filter_intervals(df, chrom, start, end)
    edges = np.linspace(start, end, n_bins + 1, dtype=int)
    mids = ((edges[:-1] + edges[1:]) / 2).astype(int)
    values = np.zeros(n_bins, dtype=float)
    if locus.empty:
        return pd.DataFrame({"x": mids, "y": values})

    for row in locus.itertuples():
        overlap_start = np.searchsorted(edges, max(start, int(row.start)), side="right") - 1
        overlap_end = np.searchsorted(edges, min(end, int(row.end)), side="left")
        overlap_start = max(0, overlap_start)
        overlap_end = min(n_bins, overlap_end)
        if overlap_end <= overlap_start:
            continue
        values[overlap_start:overlap_end] = np.maximum(values[overlap_start:overlap_end], float(row.signal))
    return pd.DataFrame({"x": mids, "y": values})


def assign_gene_rows(genes: pd.DataFrame) -> pd.DataFrame:
    if genes.empty:
        genes["track_row"] = []
        return genes
    genes = genes.sort_values(["start", "end"]).copy()
    row_ends: list[int] = []
    track_rows: list[int] = []
    for gene in genes.itertuples():
        assigned = False
        for idx, current_end in enumerate(row_ends):
            if int(gene.start) > current_end:
                row_ends[idx] = int(gene.end)
                track_rows.append(idx)
                assigned = True
                break
        if not assigned:
            row_ends.append(int(gene.end))
            track_rows.append(len(row_ends) - 1)
    genes["track_row"] = track_rows
    return genes


def classify_tads(metrics_df: pd.DataFrame) -> pd.DataFrame:
    if metrics_df.empty:
        metrics_df["tad_type"] = []
        return metrics_df

    out = metrics_df.copy()
    atac_hi = out["mean_accessibility"].quantile(0.75) if "mean_accessibility" in out else 0.0
    atac_lo = out["mean_accessibility"].quantile(0.25) if "mean_accessibility" in out else 0.0
    rna_hi = out["mean_rna"].quantile(0.75) if "mean_rna" in out else 0.0
    rna_lo = out["mean_rna"].quantile(0.25) if "mean_rna" in out else 0.0
    ctcf_hi = out["boundary_ctcf_score"].quantile(0.75) if "boundary_ctcf_score" in out else 0.0
    ctcf_lo = out["boundary_ctcf_score"].quantile(0.25) if "boundary_ctcf_score" in out else 0.0

    labels = []
    for row in out.itertuples():
        if row.mean_accessibility >= atac_hi and row.mean_rna >= rna_hi:
            labels.append("Active TAD")
        elif row.mean_accessibility <= atac_lo and row.mean_rna <= rna_lo:
            labels.append("Repressed TAD")
        elif row.boundary_ctcf_score >= ctcf_hi:
            labels.append("Boundary-strong TAD")
        elif row.boundary_ctcf_score <= ctcf_lo and row.mean_accessibility <= atac_lo:
            labels.append("Weak TAD")
        else:
            labels.append("Mixed TAD")
    out["tad_type"] = labels
    return out


def add_normalized_scores(metrics_df: pd.DataFrame) -> pd.DataFrame:
    out = metrics_df.copy()

    def _norm(col: str) -> pd.Series:
        values = pd.to_numeric(out[col], errors="coerce").fillna(0.0)
        lo = values.quantile(0.05)
        hi = values.quantile(0.95)
        if hi <= lo:
            return pd.Series(0.0, index=out.index)
        return ((values - lo) / (hi - lo)).clip(0, 1)

    out["activity_score"] = (
        0.35 * _norm("mean_accessibility") +
        0.20 * _norm("max_accessibility") +
        0.30 * _norm("mean_rna") +
        0.15 * _norm("expressed_gene_count")
    )
    out["boundary_confidence_score"] = (
        0.55 * _norm("boundary_ctcf_score") +
        0.25 * (1 - _norm("ctcf_asymmetry")) +
        0.20 * _norm("accessibility_boundary_score")
    )
    out["rna_accessibility_ratio"] = out["mean_rna"] / (out["mean_accessibility"] + 1e-9)
    return out


@st.cache_data
def compute_tad_metrics(
    tad_df: pd.DataFrame,
    boundaries_df: pd.DataFrame,
    genes_df: Optional[pd.DataFrame],
    atac_df: Optional[pd.DataFrame],
    chip_df: Optional[pd.DataFrame],
    rna_df: Optional[pd.DataFrame],
    boundary_window_bp: int = 25_000,
) -> pd.DataFrame:
    if tad_df.empty:
        return tad_df.copy()

    out = tad_df.copy().reset_index(drop=True)
    out["mean_accessibility"] = 0.0
    out["max_accessibility"] = 0.0
    out["mean_rna"] = 0.0
    out["max_rna"] = 0.0
    out["gene_count"] = 0
    out["expressed_gene_count"] = 0
    out["fraction_genes_expressed"] = 0.0
    out["gene_density_per_mb"] = 0.0
    out["accessibility_peak_count"] = 0
    out["accessibility_coverage"] = 0.0
    out["ctcf_left"] = 0.0
    out["ctcf_right"] = 0.0
    out["boundary_ctcf_score"] = 0.0
    out["ctcf_asymmetry"] = 0.0
    out["accessibility_left_boundary"] = 0.0
    out["accessibility_right_boundary"] = 0.0
    out["accessibility_boundary_score"] = 0.0
    out["top_genes"] = ""
    out["top_expressed_gene"] = ""

    for idx, row in out.iterrows():
        chrom = str(row["chrom"])
        start = int(row["start"])
        end = int(row["end"])

        if atac_df is not None:
            out.at[idx, "mean_accessibility"] = interval_weighted_mean(atac_df, chrom, start, end)
            out.at[idx, "max_accessibility"] = interval_max_signal(atac_df, chrom, start, end)
            out.at[idx, "accessibility_peak_count"] = interval_count(atac_df, chrom, start, end)
            out.at[idx, "accessibility_coverage"] = interval_coverage_fraction(atac_df, chrom, start, end)
        if rna_df is not None:
            out.at[idx, "mean_rna"] = interval_weighted_mean(rna_df, chrom, start, end)
            out.at[idx, "max_rna"] = interval_max_signal(rna_df, chrom, start, end)
        if genes_df is not None:
            genes = filter_intervals(genes_df, chrom, start, end)
            out.at[idx, "gene_count"] = int(len(genes))
            out.at[idx, "gene_density_per_mb"] = float(len(genes) / max((end - start) / 1_000_000, 1e-9))
            if not genes.empty:
                out.at[idx, "top_genes"] = ", ".join(genes["gene_name"].dropna().astype(str).head(5))
            if not genes.empty and rna_df is not None:
                genes_ranked = genes_with_rna_values(genes, rna_df, chrom)
                expressed = int((genes_ranked["rna_signal"] > 0).sum())
                out.at[idx, "expressed_gene_count"] = expressed
                out.at[idx, "fraction_genes_expressed"] = float(expressed / max(1, len(genes)))
                if not genes_ranked.empty:
                    top = genes_ranked.iloc[0]
                    out.at[idx, "top_expressed_gene"] = f"{top['gene_name']} ({float(top['rna_signal']):.2f})"

        if chip_df is not None:
            left_start = max(0, start - boundary_window_bp)
            left_end = start + boundary_window_bp
            right_start = max(0, end - boundary_window_bp)
            right_end = end + boundary_window_bp
            left_ctcf = interval_weighted_mean(chip_df, chrom, left_start, left_end)
            right_ctcf = interval_weighted_mean(chip_df, chrom, right_start, right_end)
            out.at[idx, "ctcf_left"] = left_ctcf
            out.at[idx, "ctcf_right"] = right_ctcf
            out.at[idx, "boundary_ctcf_score"] = (left_ctcf + right_ctcf) / 2.0
            out.at[idx, "ctcf_asymmetry"] = abs(left_ctcf - right_ctcf) / (left_ctcf + right_ctcf + 1e-9)

        if atac_df is not None:
            left_start = max(0, start - boundary_window_bp)
            left_end = start + boundary_window_bp
            right_start = max(0, end - boundary_window_bp)
            right_end = end + boundary_window_bp
            left_access = interval_weighted_mean(atac_df, chrom, left_start, left_end)
            right_access = interval_weighted_mean(atac_df, chrom, right_start, right_end)
            out.at[idx, "accessibility_left_boundary"] = left_access
            out.at[idx, "accessibility_right_boundary"] = right_access
            out.at[idx, "accessibility_boundary_score"] = (left_access + right_access) / 2.0

    out = classify_tads(out)
    out = add_normalized_scores(out)
    return out


def style_tad_rows(df: pd.DataFrame) -> "pd.io.formats.style.Styler":
    color_map = {
        "Active TAD": "#e6ffed",
        "Repressed TAD": "#f5f5f5",
        "Boundary-strong TAD": "#e8f1ff",
        "Weak TAD": "#fff3e0",
        "Mixed TAD": "#faf6ff",
    }

    def _row_style(row: pd.Series) -> list[str]:
        color = color_map.get(row.get("tad_type", ""), "")
        return [f"background-color: {color}"] * len(row)

    return df.style.apply(_row_style, axis=1)


def render_overview_plots(metrics_df: pd.DataFrame) -> None:
    if metrics_df.empty:
        return
    tab1, tab2, tab3 = st.tabs(["Activity", "Boundaries", "Classes"])
    with tab1:
        fig = go.Figure()
        fig.add_trace(
            go.Scatter(
                x=metrics_df["mean_accessibility"],
                y=metrics_df["mean_rna"],
                mode="markers",
                marker=dict(
                    size=np.clip(metrics_df["gene_count"] + 3, 4, 18),
                    color=metrics_df["activity_score"],
                    colorscale="Viridis",
                    colorbar=dict(title="activity"),
                    opacity=0.65,
                ),
                text=metrics_df["tad_id"],
                hovertemplate="TAD %{text}<br>access=%{x:.3f}<br>RNA=%{y:.3f}<extra></extra>",
            )
        )
        fig.update_layout(height=360, xaxis_title="Mean accessibility", yaxis_title="Mean RNA", margin=dict(l=20, r=20, t=20, b=20))
        st.plotly_chart(fig, use_container_width=True)
    with tab2:
        fig = make_subplots(rows=1, cols=2, subplot_titles=["Boundary CTCF vs Length", "CTCF Asymmetry"])
        fig.add_trace(
            go.Scatter(
                x=metrics_df["length_bp"] / 1_000_000,
                y=metrics_df["boundary_ctcf_score"],
                mode="markers",
                marker=dict(color=metrics_df["boundary_confidence_score"], colorscale="Blues", opacity=0.65),
                hovertemplate="Length=%{x:.2f} Mb<br>Boundary CTCF=%{y:.3f}<extra></extra>",
            ),
            row=1,
            col=1,
        )
        fig.add_trace(go.Histogram(x=metrics_df["ctcf_asymmetry"], nbinsx=40, marker_color="#c46a43"), row=1, col=2)
        fig.update_layout(height=360, showlegend=False, margin=dict(l=20, r=20, t=40, b=20))
        fig.update_xaxes(title_text="TAD size (Mb)", row=1, col=1)
        fig.update_yaxes(title_text="Boundary CTCF", row=1, col=1)
        fig.update_xaxes(title_text="CTCF asymmetry", row=1, col=2)
        st.plotly_chart(fig, use_container_width=True)
    with tab3:
        class_counts = metrics_df["tad_type"].value_counts().reset_index()
        class_counts.columns = ["tad_type", "count"]
        chrom_activity = metrics_df.groupby("chrom", as_index=False)["activity_score"].mean()
        fig = make_subplots(rows=1, cols=2, subplot_titles=["TAD Class Counts", "Mean Activity by Chromosome"])
        fig.add_trace(go.Bar(x=class_counts["tad_type"], y=class_counts["count"], marker_color="#4c78a8"), row=1, col=1)
        fig.add_trace(go.Bar(x=chrom_activity["chrom"], y=chrom_activity["activity_score"], marker_color="#59a14f"), row=1, col=2)
        fig.update_layout(height=380, showlegend=False, margin=dict(l=20, r=20, t=40, b=80))
        fig.update_xaxes(tickangle=35, row=1, col=1)
        fig.update_xaxes(tickangle=90, row=1, col=2)
        st.plotly_chart(fig, use_container_width=True)


def make_contact_map_trace(clr: cooler.Cooler, chrom: str, start: int, end: int) -> tuple[go.Heatmap, np.ndarray]:
    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.asarray(matrix, dtype=float)
    matrix[~np.isfinite(matrix)] = np.nan
    z = np.log10(np.clip(matrix, 1e-4, None))
    n = z.shape[0]
    x = np.linspace(start, end, n)
    trace = go.Heatmap(
        z=z,
        x=x,
        y=x,
        colorscale="RdBu",
        reversescale=True,
        colorbar=dict(title="log10 contact"),
        hovertemplate="x=%{x}<br>y=%{y}<br>log10=%{z:.2f}<extra></extra>",
    )
    return trace, matrix


def add_boundary_shapes(fig: go.Figure, row_count: int, locus_boundaries: pd.DataFrame, start: int, end: int) -> None:
    for _, boundary in locus_boundaries.iterrows():
        pos = int(boundary["start"])
        for row_idx in range(2, row_count + 1):
            fig.add_vline(
                x=pos,
                line_width=1,
                line_dash="dot",
                line_color="rgba(180,60,60,0.7)",
                row=row_idx,
                col=1,
            )
        fig.add_shape(
            type="line",
            x0=pos,
            x1=pos,
            y0=start,
            y1=end,
            line=dict(color="rgba(180,60,60,0.7)", dash="dot", width=1),
            row=1,
            col=1,
        )


def add_tad_highlight(fig: go.Figure, row_count: int, tad_start: int, tad_end: int) -> None:
    for row_idx in range(2, row_count + 1):
        fig.add_vrect(
            x0=tad_start,
            x1=tad_end,
            fillcolor="rgba(30,136,229,0.08)",
            line_width=0,
            row=row_idx,
            col=1,
        )


def render_locus_figure(
    clr: cooler.Cooler,
    chrom: str,
    locus_start: int,
    locus_end: int,
    selected_tad: pd.Series,
    locus_boundaries: pd.DataFrame,
    genes: pd.DataFrame,
    atac_track: Optional[pd.DataFrame],
    chip_track: Optional[pd.DataFrame],
    rna_track: Optional[pd.DataFrame],
) -> go.Figure:
    track_defs: list[tuple[str, str]] = [("Genes", "genes")]
    if atac_track is not None:
        track_defs.append(("ATAC", "atac"))
    if chip_track is not None:
        track_defs.append(("CTCF / ChIP", "chip"))
    if rna_track is not None:
        track_defs.append(("RNA", "rna"))

    rows = 1 + len(track_defs)
    row_heights = [0.52] + [0.12] * len(track_defs)
    fig = make_subplots(
        rows=rows,
        cols=1,
        shared_xaxes=True,
        vertical_spacing=0.03,
        row_heights=row_heights,
        subplot_titles=["Contact Map"] + [label for label, _ in track_defs],
    )

    heatmap, _ = make_contact_map_trace(clr, chrom, locus_start, locus_end)
    fig.add_trace(heatmap, row=1, col=1)

    genes_row = 2
    genes = assign_gene_rows(genes.copy()) if not genes.empty else genes
    for gene in genes.itertuples():
        y = int(gene.track_row)
        fig.add_trace(
            go.Scatter(
                x=[int(gene.start), int(gene.end)],
                y=[y, y],
                mode="lines+text",
                line=dict(width=10, color="#4c78a8"),
                text=[None, gene.gene_name],
                textposition="top center",
                hovertemplate=f"{gene.gene_name}<br>{gene.chrom}:{gene.start:,}-{gene.end:,}<extra></extra>",
                showlegend=False,
            ),
            row=genes_row,
            col=1,
        )
    if genes.empty:
        fig.add_annotation(text="No genes loaded", x=locus_start, y=0, showarrow=False, row=genes_row, col=1)

    current_row = 3
    for label, key in track_defs[1:]:
        source = {"atac": atac_track, "chip": chip_track, "rna": rna_track}[key]
        sampled = sample_signal_track(source if source is not None else pd.DataFrame(columns=["chrom", "start", "end", "signal"]), chrom, locus_start, locus_end)
        fig.add_trace(
            go.Scatter(
                x=sampled["x"],
                y=sampled["y"],
                mode="lines",
                fill="tozeroy",
                line=dict(width=1.5),
                name=label,
                hovertemplate=f"{label}<br>%{{x}}<br>%{{y:.3f}}<extra></extra>",
                showlegend=False,
            ),
            row=current_row,
            col=1,
        )
        current_row += 1

    add_tad_highlight(fig, rows, int(selected_tad["start"]), int(selected_tad["end"]))
    add_boundary_shapes(fig, rows, locus_boundaries, locus_start, locus_end)

    fig.update_layout(
        height=300 + 170 * rows,
        margin=dict(l=30, r=30, t=60, b=30),
        title=f"{chrom}:{locus_start:,}-{locus_end:,}",
    )
    fig.update_xaxes(title_text="Genomic coordinate (bp)", row=rows, col=1)
    fig.update_yaxes(scaleanchor="x", scaleratio=1, row=1, col=1, title_text="Position")
    fig.update_yaxes(showticklabels=False, row=genes_row, col=1, title_text="")
    for row_idx in range(3, rows + 1):
        fig.update_yaxes(showgrid=False, zeroline=False, row=row_idx, col=1)
    return fig


def summarize_selected_tad(
    tad_row: pd.Series,
    genes: pd.DataFrame,
    atac_track: Optional[pd.DataFrame],
    chip_track: Optional[pd.DataFrame],
    rna_track: Optional[pd.DataFrame],
) -> str:
    chrom = str(tad_row["chrom"])
    tad_start = int(tad_row["start"])
    tad_end = int(tad_row["end"])
    pieces: list[str] = []

    gene_names = genes["gene_name"].dropna().astype(str).tolist()[:5]
    pieces.append(f"This TAD spans {((tad_end - tad_start) / 1000):.0f} kb")
    if gene_names:
        pieces.append(f"and contains {len(genes)} genes ({', '.join(gene_names)})")
    else:
        pieces.append("and has no loaded genes in the visible annotation set")

    for label, track in [("ATAC", atac_track), ("CTCF/ChIP", chip_track), ("RNA", rna_track)]:
        if track is None:
            continue
        inside = filter_intervals(track, chrom, tad_start, tad_end)
        if inside.empty:
            pieces.append(f"{label} has no overlapping signal")
        else:
            pieces.append(f"{label} peaks at {inside['signal'].max():.2f} with mean interval signal {inside['signal'].mean():.2f}")
    if "tad_type" in tad_row:
        pieces.append(f"heuristic class: {tad_row['tad_type']}")
    if "boundary_ctcf_score" in tad_row:
        pieces.append(f"boundary CTCF score {float(tad_row['boundary_ctcf_score']):.2f}")
    if "expressed_gene_count" in tad_row:
        pieces.append(f"{int(tad_row['expressed_gene_count'])} expressed genes in-domain")
    return ". ".join(pieces) + "."


st.title("Genome-wide TAD Browser")
st.caption("Browse a called TAD set, then inspect a synchronized locus view with structure and optional regulatory tracks.")

manifest = load_manifest(str(DEFAULT_MANIFEST))

with st.sidebar:
    st.header("Callset")
    prefix_str = st.text_input("Callset prefix", str(DEFAULT_PREFIX))
    cooler_path = st.text_input("Cooler file", str(DEFAULT_COOLER))
    resolution = st.selectbox("Resolution", [5000, 10000, 25000, 50000, 100000, 250000], index=4)
    flank_label = st.selectbox("Locus padding", ["TAD only", "TAD ± 50 kb", "TAD ± 100 kb", "TAD ± 250 kb", "TAD ± 500 kb"], index=3)
    flank_bp = {
        "TAD only": 0,
        "TAD ± 50 kb": 50_000,
        "TAD ± 100 kb": 100_000,
        "TAD ± 250 kb": 250_000,
        "TAD ± 500 kb": 500_000,
    }[flank_label]

    st.header("Tracks")
    gene_path = st.text_input("Gene annotation path (GTF/BED)", manifest.get("genes", ""))
    atac_path = st.text_input("ATAC / accessibility path", manifest.get("atac", ""))
    chip_path = st.text_input("CTCF / ChIP path", manifest.get("ctcf", ""))
    rna_path = st.text_input("RNA path", manifest.get("rna", ""))

    st.subheader("Validation")
    for name, path_str, kind in [
        ("Genes", gene_path, "genes"),
        ("CTCF", chip_path, "signal"),
        ("ATAC", atac_path, "signal"),
        ("RNA", rna_path, "signal"),
    ]:
        if path_str and os.path.exists(path_str):
            issues = validate_track(path_str, kind=kind)
            if issues:
                st.warning(f"{name}: " + "; ".join(issues))
            else:
                st.success(f"{name}: OK")

tads, boundaries, tad_summary, prefix = load_callset(prefix_str)

genes_df: Optional[pd.DataFrame] = None
atac_df: Optional[pd.DataFrame] = None
chip_df: Optional[pd.DataFrame] = None
rna_df: Optional[pd.DataFrame] = None

track_load_errors: list[str] = []
for label, path_str, target in [
    ("Genes", gene_path, "genes"),
    ("ATAC", atac_path, "atac"),
    ("CTCF / ChIP", chip_path, "chip"),
    ("RNA", rna_path, "rna"),
]:
    if not path_str.strip():
        continue
    try:
        if target == "genes":
            genes_df = load_gene_annotations(path_str.strip())
        else:
            loaded = load_signal_intervals(path_str.strip())
            if target == "atac":
                atac_df = loaded
            elif target == "chip":
                chip_df = loaded
            else:
                rna_df = loaded
    except Exception as exc:
        track_load_errors.append(f"{label}: {exc}")

if track_load_errors:
    for err in track_load_errors:
        st.warning(err)

metrics_df = compute_tad_metrics(tad_summary, boundaries, genes_df, atac_df, chip_df, rna_df)

metric_cols = st.columns(5)
metric_cols[0].metric("TADs", f"{len(tads):,}")
metric_cols[1].metric("Chromosomes", f"{tads['chrom'].nunique():,}" if not tads.empty else "0")
metric_cols[2].metric("Median size", f"{int(tads['length_bp'].median() / 1000):,} kb" if not tads.empty else "0 kb")
metric_cols[3].metric("Strong boundaries", f"{int((boundaries['boundary_class'] == 'strong').sum()):,}" if not boundaries.empty else "0")
metric_cols[4].metric("Tracks loaded", sum(x is not None for x in [genes_df, atac_df, chip_df, rna_df]))

summary_image = prefix.with_suffix(".summary.png")
if summary_image.exists():
    st.image(str(summary_image), caption=summary_image.name, use_container_width=True)

st.subheader("Genome-wide Overview")
render_overview_plots(metrics_df)

left, right = st.columns([1.1, 0.9])

with left:
    st.subheader("TAD Table")
    chromosomes = ["All"] + sorted(tads["chrom"].unique().tolist())
    chrom_filter = st.selectbox("Chromosome", chromosomes)
    filtered = metrics_df.copy()
    if chrom_filter != "All":
        filtered = filtered[filtered["chrom"] == chrom_filter].copy()

    min_len = int(filtered["length_bp"].min())
    max_len = int(filtered["length_bp"].max())
    length_range = st.slider(
        "TAD length (bp)",
        min_value=min_len,
        max_value=max_len,
        value=(min_len, max_len),
        step=max(10_000, (max_len - min_len) // 100),
    )
    filtered = filtered[
        (filtered["length_bp"] >= length_range[0]) &
        (filtered["length_bp"] <= length_range[1])
    ].copy()
    filtered = filtered.sort_values(["chrom", "start"]).reset_index(drop=True)

    tad_type_options = ["All"] + sorted(filtered["tad_type"].dropna().unique().tolist())
    tad_type_filter = st.selectbox("TAD type", tad_type_options)
    if tad_type_filter != "All":
        filtered = filtered[filtered["tad_type"] == tad_type_filter].copy()

    sort_options = [
        "start",
        "activity_score",
        "boundary_confidence_score",
        "mean_accessibility",
        "accessibility_peak_count",
        "accessibility_coverage",
        "mean_rna",
        "boundary_ctcf_score",
        "ctcf_asymmetry",
        "gene_count",
        "expressed_gene_count",
        "fraction_genes_expressed",
        "length_bp",
    ]
    sort_by = st.selectbox("Sort by", sort_options, index=1)
    sort_desc = st.checkbox("Descending", value=sort_by != "start")
    filtered = filtered.sort_values(sort_by, ascending=not sort_desc).reset_index(drop=True)

    show_cols = [
        "chrom", "start", "end", "tad_id", "length_bp", "tad_type",
        "activity_score", "boundary_confidence_score",
        "mean_accessibility", "mean_rna", "boundary_ctcf_score", "ctcf_asymmetry",
        "accessibility_peak_count", "accessibility_coverage",
        "gene_count", "expressed_gene_count", "fraction_genes_expressed",
        "top_expressed_gene",
        "left_boundary_class", "right_boundary_class",
    ]
    show_cols += [c for c in get_assay_columns(filtered) if c.endswith("_mean")]
    show_cols = [c for c in show_cols if c in filtered.columns]
    st.dataframe(style_tad_rows(filtered[show_cols]), use_container_width=True, height=420)

    with st.expander("Interesting TADs"):
        pick_cols = [c for c in ["chrom", "start", "end", "tad_id", "tad_type", "mean_accessibility", "mean_rna", "boundary_ctcf_score", "gene_count"] if c in filtered.columns]
        if not filtered.empty:
            st.markdown("`Top active`")
            st.dataframe(filtered.sort_values(["mean_accessibility", "mean_rna"], ascending=False).head(10)[pick_cols], use_container_width=True, height=220)
            st.markdown("`Weakest boundaries`")
            st.dataframe(filtered.sort_values("boundary_ctcf_score", ascending=True).head(10)[pick_cols], use_container_width=True, height=220)
            st.markdown("`Most asymmetric boundaries`")
            st.dataframe(filtered.sort_values("ctcf_asymmetry", ascending=False).head(10)[pick_cols], use_container_width=True, height=220)
            st.markdown("`Most gene-dense`")
            st.dataframe(filtered.sort_values("gene_count", ascending=False).head(10)[pick_cols], use_container_width=True, height=220)

with right:
    st.subheader("Selected TAD")
    tad_ids = filtered["tad_id"].tolist()
    selected_tad_id = st.selectbox("Select TAD", tad_ids, index=0 if tad_ids else None)
    if selected_tad_id is not None:
        row = filtered[filtered["tad_id"] == selected_tad_id].iloc[0]
        chrom = str(row["chrom"])
        tad_start = int(row["start"])
        tad_end = int(row["end"])
        locus_start = max(0, tad_start - flank_bp)
        locus_end = tad_end + flank_bp

        clr = load_cooler(cooler_path, resolution)
        locus_boundaries = boundaries[
            (boundaries["chrom"] == chrom) &
            (boundaries["start"] >= locus_start) &
            (boundaries["end"] <= locus_end)
        ].copy()
        locus_genes = filter_intervals(genes_df, chrom, locus_start, locus_end) if genes_df is not None else pd.DataFrame(columns=["chrom", "start", "end", "gene_name", "strand"])

        meta_cols = st.columns(6)
        meta_cols[0].metric("Span", f"{(tad_end - tad_start) / 1000:.0f} kb")
        meta_cols[1].metric("Genes in view", len(locus_genes))
        meta_cols[2].metric("Boundaries in view", len(locus_boundaries))
        meta_cols[3].metric("Mean accessibility", f"{float(row.get('mean_accessibility', 0.0)):.2f}")
        meta_cols[4].metric("Boundary CTCF", f"{float(row.get('boundary_ctcf_score', 0.0)):.2f}")
        meta_cols[5].metric("Activity score", f"{float(row.get('activity_score', 0.0)):.2f}")

        st.caption(f"Centered on `{chrom}:{locus_start:,}-{locus_end:,}`")
        st.info(summarize_selected_tad(row, filter_intervals(locus_genes, chrom, tad_start, tad_end), atac_df, chip_df, rna_df))

        assay_cols = get_assay_columns(filtered)
        if assay_cols:
            assay_df = row[assay_cols].to_frame(name="value")
            assay_df.index.name = "feature"
            st.dataframe(assay_df, use_container_width=True, height=220)

        stats_cols = [c for c in [
            "tad_type", "activity_score", "boundary_confidence_score",
            "mean_accessibility", "max_accessibility", "accessibility_peak_count", "accessibility_coverage",
            "mean_rna", "max_rna", "gene_count", "expressed_gene_count", "fraction_genes_expressed",
            "top_genes", "top_expressed_gene",
            "ctcf_left", "ctcf_right", "ctcf_asymmetry", "boundary_ctcf_score",
            "accessibility_left_boundary", "accessibility_right_boundary", "accessibility_boundary_score",
        ] if c in row.index]
        if stats_cols:
            st.dataframe(row[stats_cols].to_frame(name="value"), use_container_width=True, height=250)

        genes_in_tad = filter_intervals(locus_genes, chrom, tad_start, tad_end)
        genes_in_tad = genes_with_rna_values(genes_in_tad, rna_df, chrom) if not genes_in_tad.empty else genes_in_tad
        if not genes_in_tad.empty:
            with st.expander("Genes in selected TAD", expanded=True):
                display_genes = genes_in_tad[["gene_name", "chrom", "start", "end", "strand", "rna_signal"]].head(50)
                st.dataframe(display_genes, use_container_width=True, height=260)

        with st.expander("Boundary Details", expanded=False):
            boundary_detail = pd.DataFrame(
                [
                    {
                        "side": "left",
                        "position": tad_start,
                        "ctcf": row.get("ctcf_left", 0.0),
                        "accessibility": row.get("accessibility_left_boundary", 0.0),
                        "class": row.get("left_boundary_class", ""),
                    },
                    {
                        "side": "right",
                        "position": tad_end,
                        "ctcf": row.get("ctcf_right", 0.0),
                        "accessibility": row.get("accessibility_right_boundary", 0.0),
                        "class": row.get("right_boundary_class", ""),
                    },
                ]
            )
            st.dataframe(boundary_detail, use_container_width=True)

        figure = render_locus_figure(
            clr,
            chrom,
            locus_start,
            locus_end,
            row,
            locus_boundaries,
            locus_genes,
            atac_df,
            chip_df,
            rna_df,
        )
        st.plotly_chart(figure, use_container_width=True)

st.subheader("Boundary Summary")
boundary_show = boundaries[["chrom", "start", "end", "boundary_class", "prominence"]].copy()
if chrom_filter != "All":
    boundary_show = boundary_show[boundary_show["chrom"] == chrom_filter]
st.dataframe(boundary_show, use_container_width=True, height=280)

with st.expander("What to load here"):
    st.markdown(
        """
Use this app as a stacked locus browser:

- `Gene annotation`: GTF gene rows or BED intervals
- `ATAC / CTCF / RNA`: BED, bedGraph, TSV, CSV, or bigWig if `pyBigWig` is installed
- All tracks should use the same genome assembly and chromosome naming as the cooler

The shared locus panel aligns:
- contact map
- TAD highlight
- boundary markers
- genes
- ATAC
- CTCF / ChIP
- RNA

Current bundled mm10 defaults:
- genes: GENCODE M25 on mm10
- CTCF: ENCODE E14TG2a.4 CTCF peaks
- accessibility: ENCODE ES-E14 DNase peaks used as the closest accessibility fallback
- RNA: ENCODE ES-E14 gene quantification lifted onto GENCODE gene intervals
"""
    )

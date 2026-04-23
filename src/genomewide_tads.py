from __future__ import annotations

import csv
from pathlib import Path
from typing import Optional

import cooler
import cooltools
import numpy as np
import pandas as pd
from src.tad_boundaries import call_tad_intervals, score_boundary_prominence


def _normalize_chromosome_list(
    clr: cooler.Cooler,
    chromosomes: Optional[list[str]] = None,
    min_chrom_size_bp: int = 0,
) -> list[str]:
    chromsizes = clr.chromsizes
    if chromosomes is None:
        chromosomes = list(chromsizes.index)
    return [chrom for chrom in chromosomes if chrom in chromsizes.index and int(chromsizes[chrom]) >= min_chrom_size_bp]


def compute_genomewide_insulation(
    clr: cooler.Cooler,
    window_bp: int,
    *,
    chromosomes: Optional[list[str]] = None,
    ignore_diags: int = 2,
    min_chrom_size_bp: int = 0,
) -> pd.DataFrame:
    """
    Compute a single-window insulation table across selected chromosomes.

    TAD calling still comes from the contact map; this is the genome-wide
    entry point that prepares one insulation track for downstream boundary
    and domain calling.
    """
    chroms = _normalize_chromosome_list(clr, chromosomes, min_chrom_size_bp=min_chrom_size_bp)
    if not chroms:
        raise ValueError("No chromosomes passed filtering; check chromosome names or min_chrom_size_bp.")

    view_df = pd.DataFrame(
        {
            "chrom": chroms,
            "start": [0] * len(chroms),
            "end": [int(clr.chromsizes[c]) for c in chroms],
            "name": chroms,
        }
    )

    return cooltools.insulation(
        clr,
        [window_bp],
        ignore_diags=ignore_diags,
        view_df=view_df,
    ).reset_index(drop=True)


def call_genomewide_tads(
    clr: cooler.Cooler,
    window_bp: int,
    *,
    chromosomes: Optional[list[str]] = None,
    ignore_diags: int = 2,
    prominence_thresholds: tuple[float, float] = (0.2, 0.5),
    boundary_class_min: str = "weak",
    min_chrom_size_bp: int = 0,
    min_tad_length_bp: int = 100_000,
    max_tad_length_bp: int = 3_000_000,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Call genome-wide TADs from a cooler contact map.

    Returns:
        insulation_table, boundary_table, tad_table
    """
    insulation_table = compute_genomewide_insulation(
        clr,
        window_bp,
        chromosomes=chromosomes,
        ignore_diags=ignore_diags,
        min_chrom_size_bp=min_chrom_size_bp,
    )
    boundary_table = score_boundary_prominence(
        insulation_table,
        window_bp=window_bp,
        prominence_thresholds=prominence_thresholds,
    )
    tad_table = call_tad_intervals(
        boundary_table,
        clr,
        boundary_class_min=boundary_class_min,
        max_tad_length_bp=max_tad_length_bp,
    )
    if min_tad_length_bp > 0 and not tad_table.empty:
        tad_table = tad_table[tad_table["length_bp"] >= min_tad_length_bp].copy()
        tad_table = tad_table.reset_index(drop=True)
        tad_table["tad_id"] = np.arange(len(tad_table))
    return insulation_table, boundary_table, tad_table


def load_interval_signal_table(path: str | Path) -> pd.DataFrame:
    """
    Load a BED/BEDGRAPH-like interval table and normalize to:
    chrom, start, end, signal

    Accepted conventions:
    - 4-column BEDGRAPH: chrom start end value
    - BED/narrowPeak-like tables with a signal-ish column
    - TSV/CSV with named columns
    """
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"Signal file not found: {path}")

    sep = "," if path.suffix.lower() == ".csv" else "\t"
    with path.open("r", newline="") as handle:
        sample = "".join(handle.readline() for _ in range(3))
    has_header = csv.Sniffer().has_header(sample) if sample.strip() else False
    df = pd.read_csv(path, sep=sep, comment="#", header=0 if has_header else None)
    if df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "signal"])

    if list(df.columns[:4]) == [0, 1, 2, 3]:
        df = df.rename(columns={0: "chrom", 1: "start", 2: "end", 3: "signal"})
    elif not {"chrom", "start", "end"}.issubset(df.columns):
        if df.shape[1] < 4:
            raise ValueError(f"{path} must have at least 4 columns or named chrom/start/end columns.")
        df = df.iloc[:, :4].copy()
        df.columns = ["chrom", "start", "end", "signal"]
    else:
        score_candidates = [
            "signal",
            "signal_value",
            "score",
            "value",
            "normalized_count",
            "expression",
            "tpm",
            "fpkm",
            "rpkm",
        ]
        signal_col = next((col for col in score_candidates if col in df.columns), None)
        if signal_col is None:
            extra_cols = [c for c in df.columns if c not in {"chrom", "start", "end", "name", "strand"}]
            if not extra_cols:
                raise ValueError(f"Could not find a signal column in {path}")
            signal_col = extra_cols[0]
        df = df.rename(columns={signal_col: "signal"})
        df = df[["chrom", "start", "end", "signal"]].copy()

    df["start"] = df["start"].astype(int)
    df["end"] = df["end"].astype(int)
    df["signal"] = pd.to_numeric(df["signal"], errors="coerce").fillna(0.0)
    return df.sort_values(["chrom", "start", "end"]).reset_index(drop=True)


def summarize_track_over_intervals(
    intervals_df: pd.DataFrame,
    signal_df: pd.DataFrame,
    *,
    prefix: str,
) -> pd.DataFrame:
    """
    Compute overlap-weighted signal summaries for interval features.
    """
    base = intervals_df.copy().reset_index(drop=True)
    if base.empty:
        return base

    length_bp = (base["end"] - base["start"]).clip(lower=1)
    base[f"{prefix}_mean"] = 0.0
    base[f"{prefix}_max"] = 0.0
    base[f"{prefix}_coverage"] = 0.0
    base[f"{prefix}_overlap_count"] = 0

    if signal_df.empty:
        return base

    signal_by_chrom = {
        chrom: grp[["start", "end", "signal"]].sort_values(["start", "end"]).reset_index(drop=True)
        for chrom, grp in signal_df.groupby("chrom")
    }

    for interval_id, row in base.iterrows():
        chrom = row["chrom"]
        if chrom not in signal_by_chrom:
            continue

        signals = signal_by_chrom[chrom]
        overlaps = signals[(signals["end"] > row["start"]) & (signals["start"] < row["end"])].copy()
        if overlaps.empty:
            continue

        overlap_bp = (
            np.minimum(overlaps["end"].to_numpy(), row["end"]) -
            np.maximum(overlaps["start"].to_numpy(), row["start"])
        ).clip(min=0)
        valid = overlap_bp > 0
        if not np.any(valid):
            continue

        overlap_bp = overlap_bp[valid]
        signals_overlap = overlaps.loc[valid, "signal"].to_numpy(dtype=float)
        base.loc[interval_id, f"{prefix}_mean"] = float(np.sum(overlap_bp * signals_overlap) / length_bp.iloc[interval_id])
        base.loc[interval_id, f"{prefix}_max"] = float(np.max(signals_overlap))
        base.loc[interval_id, f"{prefix}_coverage"] = float(np.sum(overlap_bp) / length_bp.iloc[interval_id])
        base.loc[interval_id, f"{prefix}_overlap_count"] = int(np.sum(valid))

    return base


def make_boundary_windows(
    tad_df: pd.DataFrame,
    *,
    chromsizes: pd.Series,
    flank_bp: int = 25_000,
    include_chrom_edges: bool = False,
) -> pd.DataFrame:
    """
    Build boundary-centered windows from TAD intervals.
    """
    if tad_df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "boundary_bp", "boundary_type"])

    rows: list[dict[str, object]] = []
    for _, row in tad_df.iterrows():
        chrom = row["chrom"]
        chrom_len = int(chromsizes[chrom])
        left = int(row["start"])
        right = int(row["end"])

        if include_chrom_edges or left > 0:
            rows.append({"chrom": chrom, "boundary_bp": left, "boundary_type": "left"})
        if include_chrom_edges or right < chrom_len:
            rows.append({"chrom": chrom, "boundary_bp": right, "boundary_type": "right"})

    boundary_df = pd.DataFrame(rows).drop_duplicates(["chrom", "boundary_bp"]).reset_index(drop=True)
    if boundary_df.empty:
        return pd.DataFrame(columns=["chrom", "start", "end", "boundary_bp", "boundary_type"])

    boundary_df["start"] = (boundary_df["boundary_bp"] - flank_bp).clip(lower=0)
    boundary_df["end"] = [
        min(int(chromsizes[row.chrom]), int(row.boundary_bp + flank_bp))
        for row in boundary_df.itertuples()
    ]
    return boundary_df[["chrom", "start", "end", "boundary_bp", "boundary_type"]]


def annotate_tads_with_assays(
    tad_df: pd.DataFrame,
    *,
    chromsizes: pd.Series,
    atac_df: Optional[pd.DataFrame] = None,
    chip_df: Optional[pd.DataFrame] = None,
    rna_df: Optional[pd.DataFrame] = None,
    boundary_window_bp: int = 25_000,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Annotate called TADs and their boundaries with auxiliary assays.

    Assays are used as annotations of the called domains, not as the
    source of the TAD calls themselves.
    """
    tad_summary = tad_df.copy().reset_index(drop=True)
    boundary_summary = make_boundary_windows(
        tad_summary,
        chromsizes=chromsizes,
        flank_bp=boundary_window_bp,
        include_chrom_edges=False,
    )

    assay_tables = [
        ("atac", atac_df),
        ("chip", chip_df),
        ("rna", rna_df),
    ]
    for prefix, table in assay_tables:
        if table is None:
            table = pd.DataFrame(columns=["chrom", "start", "end", "signal"])
        tad_summary = summarize_track_over_intervals(tad_summary, table, prefix=prefix)
        boundary_summary = summarize_track_over_intervals(boundary_summary, table, prefix=f"{prefix}_boundary")

    return tad_summary, boundary_summary

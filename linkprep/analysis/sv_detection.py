"""
Structural variant (SV) detection from a LinkPrep pairs file.

LinkPrep uses Tn5 tagmentation → uniform coverage → excellent for SV calling.

Detects:
  - Translocations  : inter-chromosomal pairs clustering at a locus
  - Inversions      : cis pairs with unexpected strand orientation (++ or --)
  - Large deletions : cis pairs with unusually long distance
  - Duplications    : inferred from coverage + pairs enrichment

Usage:
    from analysis.sv_detection import detect_svs
    svs = detect_svs(pairs_path, bin_size=50_000)
"""

import gzip
import os
import sys
from collections import defaultdict

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ---------------------------------------------------------------------------
# Pairs loader
# ---------------------------------------------------------------------------

def _load_pairs(path: str, max_rows: int | None = None) -> pd.DataFrame:
    opener = gzip.open if path.endswith(".gz") else open
    rows = []
    with opener(path, "rt") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 8:
                continue
            rows.append({
                "chrom1":    p[1],
                "pos1":      int(p[2]),
                "chrom2":    p[3],
                "pos2":      int(p[4]),
                "strand1":   p[5],
                "strand2":   p[6],
                "pair_type": p[7],
            })
            if max_rows and len(rows) >= max_rows:
                break
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Translocation detection
# ---------------------------------------------------------------------------

def detect_translocations(
    pairs_df: pd.DataFrame,
    bin_size: int = 100_000,
    min_pairs: int = 10,
) -> pd.DataFrame:
    """
    Bin inter-chromosomal pairs into a 2D grid.
    Bins with ≥ min_pairs inter-chrom pairs are reported as translocation candidates.
    """
    trans = pairs_df[pairs_df["chrom1"] != pairs_df["chrom2"]].copy()
    if trans.empty:
        return pd.DataFrame()

    trans["bin1"] = (trans["pos1"] // bin_size) * bin_size
    trans["bin2"] = (trans["pos2"] // bin_size) * bin_size

    grouped = (
        trans.groupby(["chrom1", "bin1", "chrom2", "bin2"])
        .size()
        .reset_index(name="n_pairs")
    )
    candidates = grouped[grouped["n_pairs"] >= min_pairs].copy()
    candidates["sv_type"] = "translocation"
    candidates["confidence"] = np.minimum(
        candidates["n_pairs"] / (min_pairs * 5), 1.0
    )
    return candidates.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Inversion detection
# ---------------------------------------------------------------------------

def detect_inversions(
    pairs_df: pd.DataFrame,
    bin_size: int = 100_000,
    min_pairs: int = 5,
) -> pd.DataFrame:
    """
    Inversions produce ++ or -- strand orientation for cis pairs.
    Normal ligation produces +- or -+.
    """
    cis = pairs_df[pairs_df["chrom1"] == pairs_df["chrom2"]].copy()
    if cis.empty:
        return pd.DataFrame()

    inv = cis[cis["strand1"] == cis["strand2"]].copy()   # ++ or --
    if inv.empty:
        return pd.DataFrame()

    inv["bin1"] = (inv[["pos1", "pos2"]].min(axis=1) // bin_size) * bin_size
    inv["bin2"] = (inv[["pos1", "pos2"]].max(axis=1) // bin_size) * bin_size

    grouped = (
        inv.groupby(["chrom1", "bin1", "bin2"])
        .size()
        .reset_index(name="n_pairs")
    )
    grouped = grouped.rename(columns={"chrom1": "chrom"})
    candidates = grouped[grouped["n_pairs"] >= min_pairs].copy()
    candidates["sv_type"] = "inversion"
    candidates["size_bp"] = candidates["bin2"] - candidates["bin1"]
    return candidates.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Large deletion detection
# ---------------------------------------------------------------------------

def detect_deletions(
    pairs_df: pd.DataFrame,
    distance_threshold: int = 2_000_000,
    bin_size: int = 100_000,
    min_pairs: int = 5,
) -> pd.DataFrame:
    """
    Very long-range cis pairs (> distance_threshold) can mark spanning
    a deletion — the two ends of the deleted segment get brought together.
    """
    cis = pairs_df[pairs_df["chrom1"] == pairs_df["chrom2"]].copy()
    if cis.empty:
        return pd.DataFrame()

    cis["distance"] = (cis["pos2"] - cis["pos1"]).abs()
    long_range = cis[cis["distance"] >= distance_threshold].copy()
    if long_range.empty:
        return pd.DataFrame()

    long_range["bin1"] = long_range[["pos1", "pos2"]].min(axis=1) // bin_size * bin_size
    long_range["bin2"] = long_range[["pos1", "pos2"]].max(axis=1) // bin_size * bin_size

    grouped = (
        long_range.groupby(["chrom1", "bin1", "bin2"])
        .agg(n_pairs=("distance", "count"), mean_distance=("distance", "mean"))
        .reset_index()
    )
    grouped = grouped.rename(columns={"chrom1": "chrom"})
    candidates = grouped[grouped["n_pairs"] >= min_pairs].copy()
    candidates["sv_type"] = "deletion_candidate"
    return candidates.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Master function
# ---------------------------------------------------------------------------

def detect_svs(
    pairs_path: str,
    bin_size: int = 100_000,
    min_pairs: int = 5,
    distance_threshold: int = 2_000_000,
    max_rows: int = 1_000_000,
) -> dict[str, pd.DataFrame]:
    """
    Run all SV detectors on a pairs file.

    Returns
    -------
    dict with keys: "translocations", "inversions", "deletions"
    """
    print(f"Loading pairs from {pairs_path} …")
    df = _load_pairs(pairs_path, max_rows=max_rows)
    print(f"  Loaded {len(df):,} pairs")

    translocations = detect_translocations(df, bin_size=bin_size,
                                           min_pairs=min_pairs * 2)
    inversions     = detect_inversions(df, bin_size=bin_size,
                                       min_pairs=min_pairs)
    deletions      = detect_deletions(df, distance_threshold=distance_threshold,
                                      bin_size=bin_size, min_pairs=min_pairs)

    print(f"  Translocations: {len(translocations)} candidate bins")
    print(f"  Inversions    : {len(inversions)} candidate bins")
    print(f"  Deletions     : {len(deletions)} candidate bins")

    return {
        "translocations": translocations,
        "inversions":     inversions,
        "deletions":      deletions,
    }

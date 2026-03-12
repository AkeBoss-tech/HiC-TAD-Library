"""
Capture Hi-C loop calling using a CHiCAGO-inspired scoring model.

Reference: Dovetail Analysis — Capture Hi-C Loop Calling
https://dovetail-analysis.readthedocs.io/en/latest/capture/loop_calling.html

CHiCAGO models background signal as a combination of:
  - Brownian noise: negative-binomial distributed, distance-decaying signal
  - Technical noise: Poisson-distributed, distance-independent

An interaction is called significant when its score (−log weighted p-value) > 5.

This module implements:
  1. load_baitmap(path)          — parse .baitmap file (bait fragment definitions)
  2. load_rmap(path)             — parse .rmap file (all restriction/bin fragments)
  3. load_ibed(path)             — parse CHiCAGO .ibed output
  4. chicago_score(counts, bg)   — simplified CHiCAGO scoring
  5. call_loops(contacts_df)     — full loop-calling pipeline on a contact DataFrame
  6. filter_loops(loops_df)      — apply standard post-processing filters
  7. generate_sample_capture_data() — create synthetic bait + interaction data

File formats
------------
.baitmap  : chrom  start  end  bait_id  bait_name  (tab-separated)
.rmap     : chrom  start  end  frag_id              (tab-separated)
.ibed     : bait_chr  bait_start  bait_end  bait_name
            oe_chr  oe_start  oe_end  oe_name
            n_reads  score                           (10 columns)
"""

import os
import sys
from typing import Optional

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ---------------------------------------------------------------------------
# File parsers
# ---------------------------------------------------------------------------

def load_baitmap(path: str) -> pd.DataFrame:
    """Parse a .baitmap file into a DataFrame."""
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "bait_id", "bait_name"],
    )
    return df


def load_rmap(path: str) -> pd.DataFrame:
    """Parse a .rmap file into a DataFrame."""
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["chrom", "start", "end", "frag_id"],
    )
    return df


def load_ibed(path: str) -> pd.DataFrame:
    """
    Parse a CHiCAGO .ibed file.

    Columns: bait_chr, bait_start, bait_end, bait_name,
             oe_chr, oe_start, oe_end, oe_name, n_reads, score
    """
    df = pd.read_csv(
        path, sep="\t", header=None,
        names=["bait_chr", "bait_start", "bait_end", "bait_name",
               "oe_chr", "oe_start", "oe_end", "oe_name",
               "n_reads", "score"],
    )
    return df


# ---------------------------------------------------------------------------
# CHiCAGO-inspired scoring
# ---------------------------------------------------------------------------

def _brownian_background(distance_bp: np.ndarray,
                         alpha: float = 0.5,
                         beta: float = 1e-5) -> np.ndarray:
    """
    Distance-decaying Brownian noise estimate.

    bg_brownian(d) = alpha * (beta * d)^(-1)

    In the full CHiCAGO model this is fitted from data; here we use a
    sensible default that matches the typical decay slope in Micro-C / Omni-C.
    """
    d = np.maximum(distance_bp, 1)
    return alpha / (beta * d)


def _technical_noise(n_total_pairs: int, n_bins: int) -> float:
    """Flat Poisson technical noise per bin."""
    return n_total_pairs / n_bins * 0.01   # ~1% technical background


def chicago_score(
    counts: np.ndarray,
    distances: np.ndarray,
    n_total_pairs: int = 1_000_000,
    n_bins: int = 10_000,
) -> np.ndarray:
    """
    Simplified CHiCAGO score: −log10(combined p-value).

    Combines:
      - Brownian background (negative-binomial, distance-decaying)
      - Technical noise (Poisson, flat)

    Parameters
    ----------
    counts    : observed read counts per interaction
    distances : genomic distances (bp) for each interaction
    n_total_pairs, n_bins : library size parameters for noise estimation

    Returns
    -------
    scores : array of CHiCAGO-like scores (higher = more significant)
    """
    brownian  = _brownian_background(distances)
    technical = _technical_noise(n_total_pairs, n_bins)
    expected  = brownian + technical

    # Negative-binomial p-value for each interaction
    # Fit NB: mean = expected, dispersion = 0.5 (typical value from CHiCAGO paper)
    dispersion = 0.5
    r = 1.0 / dispersion   # NB size parameter
    p = r / (r + expected)

    # P(X >= counts) = 1 - CDF(counts - 1)
    p_values = 1.0 - stats.nbinom.cdf(counts - 1, r, p)
    p_values = np.clip(p_values, 1e-300, 1.0)

    scores = -np.log10(p_values)
    return scores


# ---------------------------------------------------------------------------
# Loop calling pipeline
# ---------------------------------------------------------------------------

def call_loops(
    contacts_df: pd.DataFrame,
    cutoff: float = 5.0,
    n_total_pairs: int = 1_000_000,
) -> pd.DataFrame:
    """
    Call significant capture Hi-C loops from a contacts DataFrame.

    Input DataFrame must have columns:
      bait_chr, bait_start, bait_end, bait_name,
      oe_chr, oe_start, oe_end, oe_name, n_reads

    Returns DataFrame with an added 'score' column; rows passing the
    CHiCAGO cutoff (score ≥ cutoff) are returned.
    """
    df = contacts_df.copy()

    bait_mid = (df["bait_start"] + df["bait_end"]) / 2
    oe_mid   = (df["oe_start"]   + df["oe_end"])   / 2

    # Only compute scores for cis interactions
    cis_mask = df["bait_chr"] == df["oe_chr"]
    distances = np.where(cis_mask, (oe_mid - bait_mid).abs(), np.inf)

    n_bins = len(df)
    df["distance_bp"] = distances.astype(int)
    df["score"] = chicago_score(
        df["n_reads"].values.astype(float),
        distances,
        n_total_pairs=n_total_pairs,
        n_bins=n_bins,
    )

    # Significance tiers (mirrors CHiCAGO colour scheme)
    df["significance"] = pd.cut(
        df["score"],
        bins=[-np.inf, 3, cutoff, np.inf],
        labels=["weak", "moderate", "significant"],
    )

    return df[df["score"] >= 3].reset_index(drop=True)


def filter_loops(loops_df: pd.DataFrame,
                 min_dist: int = 10_000,
                 max_dist: int = 2_000_000,
                 min_reads: int = 5,
                 cis_only: bool = True) -> pd.DataFrame:
    """
    Standard Dovetail post-processing filters:
      - Remove trans interactions
      - Remove interactions < 10 kb or > 2 Mb
      - Remove interactions with < 5 supporting reads
    """
    df = loops_df.copy()
    if cis_only:
        df = df[df["bait_chr"] == df["oe_chr"]]
    df = df[df["distance_bp"] >= min_dist]
    df = df[df["distance_bp"] <= max_dist]
    df = df[df["n_reads"]     >= min_reads]
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Sample data generator
# ---------------------------------------------------------------------------

def generate_sample_capture_data(
    rng: np.random.Generator | None = None,
    n_baits: int = 20,
    n_oe_per_bait: int = 50,
    bin_size: int = 10_000,
    chrom_sizes: dict | None = None,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Generate synthetic baitmap, rmap, and contacts DataFrames.

    Returns (baitmap_df, rmap_df, contacts_df)
    """
    if rng is None:
        rng = np.random.default_rng(42)
    if chrom_sizes is None:
        chrom_sizes = {"chr12": 120_000_000, "chr13": 100_000_000}

    # --- baitmap: evenly spaced baits per chromosome ---
    bait_rows = []
    bait_id   = 0
    for chrom, length in chrom_sizes.items():
        positions = np.linspace(5_000_000, length - 5_000_000,
                                n_baits // len(chrom_sizes), dtype=int)
        for pos in positions:
            bait_rows.append({
                "chrom":     chrom,
                "start":     int(pos),
                "end":       int(pos) + 2000,
                "bait_id":   bait_id,
                "bait_name": f"bait_{bait_id}",
            })
            bait_id += 1
    baitmap_df = pd.DataFrame(bait_rows)

    # --- rmap: 10 kb bins across genome ---
    rmap_rows = []
    frag_id   = 0
    for chrom, length in chrom_sizes.items():
        starts = np.arange(0, length, bin_size)
        for s in starts:
            rmap_rows.append({
                "chrom":   chrom,
                "start":   int(s),
                "end":     int(min(s + bin_size, length)),
                "frag_id": frag_id,
            })
            frag_id += 1
    rmap_df = pd.DataFrame(rmap_rows)

    # --- contacts: each bait × nearby OEF bins ---
    contact_rows = []
    for _, bait in baitmap_df.iterrows():
        chrom  = bait["chrom"]
        length = chrom_sizes[chrom]
        bait_mid = (bait["start"] + bait["end"]) // 2

        # Sample OEF positions: mostly nearby (distance-decay), a few distal
        near_dist  = np.abs(rng.normal(0, 200_000, n_oe_per_bait * 3)).astype(int)
        oe_offsets = np.concatenate([near_dist,
                                     rng.integers(500_000, 1_500_000, 5)])
        oe_offsets = oe_offsets[:n_oe_per_bait]
        signs      = rng.choice([-1, 1], size=n_oe_per_bait)
        oe_mids    = np.clip(bait_mid + signs * oe_offsets, 0, length - bin_size)

        for oe_mid in oe_mids:
            oe_start = int((oe_mid // bin_size) * bin_size)
            oe_end   = oe_start + bin_size
            dist     = abs(oe_mid - bait_mid)

            # Poisson counts with distance decay
            lam = max(0.5, 30 * np.exp(-dist / 300_000))
            n_reads = rng.poisson(lam)

            contact_rows.append({
                "bait_chr":   chrom,
                "bait_start": int(bait["start"]),
                "bait_end":   int(bait["end"]),
                "bait_name":  bait["bait_name"],
                "oe_chr":     chrom,
                "oe_start":   oe_start,
                "oe_end":     oe_end,
                "oe_name":    f"oe_{oe_start}",
                "n_reads":    int(n_reads),
            })

    # Inject a handful of strong loop interactions
    for _, bait in baitmap_df.iterrows():
        chrom = bait["chrom"]
        length = chrom_sizes[chrom]
        bait_mid = (bait["start"] + bait["end"]) // 2
        # Strong loop at ~100–300 kb away
        loop_dist = rng.integers(100_000, 300_000)
        oe_start  = int(np.clip(
            ((bait_mid + loop_dist) // bin_size) * bin_size,
            0, length - bin_size,
        ))
        contact_rows.append({
            "bait_chr":   chrom,
            "bait_start": int(bait["start"]),
            "bait_end":   int(bait["end"]),
            "bait_name":  bait["bait_name"],
            "oe_chr":     chrom,
            "oe_start":   oe_start,
            "oe_end":     oe_start + bin_size,
            "oe_name":    f"oe_{oe_start}_loop",
            "n_reads":    int(rng.integers(20, 80)),
        })

    contacts_df = pd.DataFrame(contact_rows)
    return baitmap_df, rmap_df, contacts_df

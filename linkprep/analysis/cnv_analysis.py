"""
Copy-number variation (CNV) analysis from per-bin coverage.

LinkPrep's uniform Tn5 fragmentation makes it suitable for CNV calling
without GC-correction artifacts typical of restriction-enzyme Hi-C.

Functions:
  load_coverage(bed_path)               — read coverage BED into DataFrame
  normalize_coverage(df)                — GC+mappability normalisation stub
  call_cnv_segments(df, ...)            — circular binary segmentation (CBS)
  assign_copy_numbers(segments, ploidy) — convert log2-ratio to copy number

Usage:
    from analysis.cnv_analysis import load_coverage, call_cnv_segments
    cov = load_coverage("coverage.bed")
    segs = call_cnv_segments(cov)
"""

import os
import sys

import numpy as np
import pandas as pd
from scipy.ndimage import uniform_filter1d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ---------------------------------------------------------------------------
# Load coverage
# ---------------------------------------------------------------------------

def load_coverage(bed_path: str) -> pd.DataFrame:
    """
    Read a 4-column BED file: chrom  start  end  coverage
    Returns DataFrame with those four columns plus log2_ratio.
    """
    df = pd.read_csv(
        bed_path,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "coverage"],
    )
    # Avoid log(0)
    df["coverage"] = df["coverage"].clip(lower=0.01)
    genome_median = df["coverage"].median()
    df["log2_ratio"] = np.log2(df["coverage"] / genome_median)
    return df


# ---------------------------------------------------------------------------
# Normalise (stub — add GC/mappability correction here)
# ---------------------------------------------------------------------------

def normalize_coverage(df: pd.DataFrame) -> pd.DataFrame:
    """
    Placeholder for GC-bias and mappability correction.

    With real data you would:
      1. Compute GC content per bin from reference FASTA
      2. Fit LOESS(coverage ~ GC) and subtract trend
      3. Mask low-mappability bins (e.g. ENCODE blacklist)

    For now we just median-normalise per chromosome.
    """
    df = df.copy()
    for chrom, grp in df.groupby("chrom"):
        med = grp["log2_ratio"].median()
        df.loc[grp.index, "log2_ratio"] -= med
    return df


# ---------------------------------------------------------------------------
# Circular Binary Segmentation (lightweight implementation)
# ---------------------------------------------------------------------------

def _cbs_segment(values: np.ndarray, min_segment: int = 5) -> list[tuple[int, int, float]]:
    """
    Simplified CBS: recursive binary split on maximum t-statistic.
    Returns list of (start_idx, end_idx, mean_value) segments.
    """
    from scipy.stats import ttest_ind

    def _split(arr, lo, hi):
        if hi - lo < min_segment * 2:
            return [(lo, hi, float(np.mean(arr[lo:hi])))]

        best_t = 0.0
        best_k = lo + min_segment
        for k in range(lo + min_segment, hi - min_segment):
            left  = arr[lo:k]
            right = arr[k:hi]
            if len(left) < 2 or len(right) < 2:
                continue
            t, p = ttest_ind(left, right)
            if abs(t) > best_t:
                best_t = abs(t)
                best_k = k

        # Stop splitting if difference is not meaningful
        if best_t < 2.0 or best_k == lo + min_segment:
            return [(lo, hi, float(np.mean(arr[lo:hi])))]

        return _split(arr, lo, best_k) + _split(arr, best_k, hi)

    return _split(values, 0, len(values))


def call_cnv_segments(
    df: pd.DataFrame,
    smoothing_bins: int = 5,
    min_segment_bins: int = 5,
) -> pd.DataFrame:
    """
    Call CNV segments per chromosome.

    Parameters
    ----------
    df               : output of load_coverage() / normalize_coverage()
    smoothing_bins   : rolling mean window before segmentation
    min_segment_bins : minimum bins per segment

    Returns
    -------
    DataFrame: chrom, start, end, mean_log2, n_bins, cn_state
    """
    segments = []
    for chrom, grp in df.groupby("chrom"):
        grp = grp.reset_index(drop=True)
        vals = grp["log2_ratio"].values.astype(float)
        # Smooth
        smoothed = uniform_filter1d(vals, size=smoothing_bins, mode="reflect")
        raw_segs = _cbs_segment(smoothed, min_segment=min_segment_bins)

        for lo, hi, mean_val in raw_segs:
            segments.append({
                "chrom":      chrom,
                "start":      int(grp["start"].iloc[lo]),
                "end":        int(grp["end"].iloc[hi - 1]),
                "mean_log2":  round(mean_val, 4),
                "n_bins":     hi - lo,
            })

    segs_df = pd.DataFrame(segments)
    if segs_df.empty:
        return segs_df

    segs_df["cn_state"] = _infer_cn_state(segs_df["mean_log2"].values)
    return segs_df


def _infer_cn_state(log2_ratios: np.ndarray, ploidy: int = 2) -> list[str]:
    """Map log2 ratio to CN state label."""
    states = []
    for lr in log2_ratios:
        if lr < -0.8:
            states.append("deletion")
        elif lr < -0.3:
            states.append("loss")
        elif lr > 0.4:
            states.append("gain")
        elif lr > 0.8:
            states.append("amplification")
        else:
            states.append("neutral")
    return states


# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------

def summarise_cnv(segments: pd.DataFrame) -> pd.DataFrame:
    """Return only non-neutral segments with size annotation."""
    non_neutral = segments[segments["cn_state"] != "neutral"].copy()
    non_neutral["size_mb"] = (non_neutral["end"] - non_neutral["start"]) / 1e6
    return non_neutral.reset_index(drop=True)

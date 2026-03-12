"""
Contact-matrix analysis for LinkPrep data.

Functions:
  load_matrix(clr, region)          — extract balanced matrix for a region
  compute_insulation(clr, region)   — insulation scores at multiple scales
  call_tads(clr, region)            — TAD intervals from insulation minima
  compute_compartments(clr, chrom)  — A/B compartments via eigenvector
  call_loops(clr, region)           — chromatin loops (donut enrichment)

All functions take a cooler.Cooler and a coordinate string
"chr12:26,000,000-28,000,000".
"""

import os
import sys

import cooler
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_region(region: str):
    chrom, rest = region.split(":")
    start_str, end_str = rest.split("-")
    return chrom, int(start_str.replace(",", "")), int(end_str.replace(",", ""))


def load_matrix(clr: cooler.Cooler, region: str) -> np.ndarray:
    """Return the balanced (ICE-normalised) contact matrix for a region.
    Falls back to raw counts if balancing weights are absent."""
    chrom, start, end = _parse_region(region)
    has_weights = "weight" in clr.bins().columns
    try:
        mat = clr.matrix(balance=has_weights).fetch(f"{chrom}:{start}-{end}")
    except ValueError:
        mat = clr.matrix(balance=False).fetch(f"{chrom}:{start}-{end}")
    return np.nan_to_num(mat.astype(float), nan=0.0)


# ---------------------------------------------------------------------------
# Insulation score & TAD calling
# ---------------------------------------------------------------------------

def compute_insulation(
    clr: cooler.Cooler,
    region: str,
    windows_bp: list | None = None,
) -> pd.DataFrame:
    """
    Sliding-diamond insulation score (Crane et al. 2015).

    Returns a DataFrame with columns:
      chrom, start, end, insulation_<w> for each window,
      boundary_strength_<w>, is_boundary_<w>
    """
    if windows_bp is None:
        windows_bp = [50_000, 100_000, 200_000]

    chrom, start, end = _parse_region(region)
    res = clr.binsize
    mat = load_matrix(clr, region)
    n   = mat.shape[0]

    bins = clr.bins().fetch(f"{chrom}:{start}-{end}").reset_index(drop=True)
    result = bins[["chrom", "start", "end"]].copy()

    for w_bp in windows_bp:
        w = max(1, w_bp // res)
        insulation = np.full(n, np.nan)
        for i in range(w, n - w):
            block = mat[i - w:i, i:i + w]
            insulation[i] = np.nanmean(block) if block.size else np.nan

        # Log-transform and normalise
        ins = np.log2(insulation + 1e-9)
        ins_mean = np.nanmean(ins)
        ins = ins - ins_mean

        # Boundary strength = prominence of the local minimum
        neg_ins = -ins
        neg_ins = np.nan_to_num(neg_ins, nan=-1e9)
        peaks, props = find_peaks(neg_ins, prominence=0.1)
        strength = np.zeros(n)
        if len(peaks):
            strength[peaks] = props["prominences"]

        is_boundary = strength > np.nanpercentile(strength[strength > 0], 50) \
            if (strength > 0).any() else np.zeros(n, bool)

        col = f"w{w_bp // 1000}k"
        result[f"insulation_{col}"]  = ins
        result[f"boundary_strength_{col}"] = strength
        result[f"is_boundary_{col}"] = is_boundary

    return result


def call_tads(insulation_df: pd.DataFrame, window_key: str | None = None) -> pd.DataFrame:
    """
    Convert boundary positions into TAD intervals.

    Parameters
    ----------
    insulation_df : output of compute_insulation()
    window_key    : e.g. "w100k"; if None, uses first available window

    Returns DataFrame with columns: chrom, start, end, tad_id
    """
    if window_key is None:
        cols = [c for c in insulation_df.columns if c.startswith("is_boundary_")]
        if not cols:
            raise ValueError("No boundary columns in insulation_df")
        window_key = cols[0].replace("is_boundary_", "")

    boundary_col = f"is_boundary_{window_key}"
    boundary_mask = insulation_df[boundary_col].values.astype(bool)

    tads = []
    chrom = insulation_df["chrom"].iloc[0]
    boundary_idx = np.where(boundary_mask)[0]
    starts = np.concatenate([[0], boundary_idx + 1])
    ends   = np.concatenate([boundary_idx, [len(insulation_df) - 1]])

    for tad_id, (s, e) in enumerate(zip(starts, ends)):
        if s >= e:
            continue
        tads.append({
            "chrom":  chrom,
            "start":  insulation_df["start"].iloc[s],
            "end":    insulation_df["end"].iloc[e],
            "tad_id": tad_id,
            "size_kb": (insulation_df["end"].iloc[e]
                        - insulation_df["start"].iloc[s]) / 1000,
        })

    return pd.DataFrame(tads)


# ---------------------------------------------------------------------------
# A/B Compartments
# ---------------------------------------------------------------------------

def compute_compartments(
    clr: cooler.Cooler,
    region: str,
) -> pd.DataFrame:
    """
    Compute A/B compartments via the first eigenvector of the
    observed/expected (OE) contact matrix.

    Returns DataFrame: chrom, start, end, E1, compartment (A/B)
    """
    chrom, start, end = _parse_region(region)
    mat = load_matrix(clr, region)
    n   = mat.shape[0]

    # Observed / expected: divide by mean contact at each diagonal
    oe = mat.copy()
    for d in range(n):
        diag = np.diag(mat, d)
        mean_d = np.nanmean(diag)
        if mean_d > 0:
            new_diag = diag / mean_d
            idx = np.arange(n - d)
            oe[idx, idx + d] = new_diag
            if d > 0:
                oe[idx + d, idx] = new_diag

    # Pearson correlation matrix
    valid = np.where(np.nansum(oe, axis=1) > 0)[0]
    corr = np.full((n, n), np.nan)
    if len(valid) > 1:
        sub = oe[np.ix_(valid, valid)]
        sub = np.nan_to_num(sub, nan=0.0)
        corr_sub = np.corrcoef(sub)
        corr[np.ix_(valid, valid)] = corr_sub

    # First eigenvector
    corr_filled = np.nan_to_num(corr, nan=0.0)
    try:
        eigenvalues, eigenvectors = np.linalg.eigh(corr_filled)
        e1 = eigenvectors[:, -1]   # largest eigenvalue
    except np.linalg.LinAlgError:
        e1 = np.zeros(n)

    # Sign convention: positive E1 → active (A) compartment (higher mean)
    # Use mean matrix row-sum as proxy for gene density
    row_sums = mat.sum(axis=1)
    if e1[row_sums > np.median(row_sums)].mean() < 0:
        e1 = -e1

    bins = clr.bins().fetch(f"{chrom}:{start}-{end}").reset_index(drop=True)
    result = bins[["chrom", "start", "end"]].copy()
    result["E1"] = e1
    result["compartment"] = np.where(e1 >= 0, "A", "B")
    return result


# ---------------------------------------------------------------------------
# Loop calling (simple donut enrichment)
# ---------------------------------------------------------------------------

def call_loops(
    clr: cooler.Cooler,
    region: str,
    inner_r: int = 3,
    outer_r: int = 7,
    min_enrichment: float = 1.5,
    fdr_threshold: float = 0.1,
) -> pd.DataFrame:
    """
    Identify chromatin loops using a donut background model.

    Returns DataFrame: chrom, start1, end1, start2, end2, enrichment, pvalue
    """
    from scipy import ndimage, stats as scipy_stats

    chrom, start, end = _parse_region(region)
    mat = load_matrix(clr, region)
    n   = mat.shape[0]
    res = clr.binsize

    # Build donut kernel
    size = 2 * outer_r + 1
    kernel = np.zeros((size, size), dtype=float)
    c = outer_r
    for i in range(size):
        for j in range(size):
            d = max(abs(i - c), abs(j - c))
            if inner_r < d <= outer_r:
                kernel[i, j] = 1.0
    kernel /= kernel.sum() or 1.0

    bg = ndimage.convolve(mat, kernel, mode="constant", cval=0.0)
    bg = np.maximum(bg, 1e-9)
    enrichment = mat / bg

    # Z-score per diagonal
    enrichment_ut = np.triu(enrichment, k=inner_r + 1)
    rows, cols = np.triu_indices(n, k=inner_r + 1)
    enrichments = enrichment_ut[rows, cols]

    z = (enrichments - enrichments.mean()) / (enrichments.std() + 1e-9)
    pvals = scipy_stats.norm.sf(z)

    # BH correction
    order = np.argsort(pvals)
    m = len(pvals)
    pvals_sorted = pvals[order]
    bh_thresholds = np.arange(1, m + 1) / m * fdr_threshold
    significant = pvals_sorted <= bh_thresholds
    if significant.any():
        max_sig = np.where(significant)[0].max()
        sig_mask_sorted = np.zeros(m, bool)
        sig_mask_sorted[:max_sig + 1] = True
    else:
        sig_mask_sorted = np.zeros(m, bool)

    sig_mask = np.zeros(m, bool)
    sig_mask[order] = sig_mask_sorted

    keep = sig_mask & ( enrichments >= min_enrichment)

    bins = clr.bins().fetch(f"{chrom}:{start}-{end}").reset_index(drop=True)
    loops = []
    for i, j, enr, pv in zip(rows[keep], cols[keep],
                                enrichments[keep], pvals[keep]):
        loops.append({
            "chrom":       chrom,
            "start1":      bins["start"].iloc[i],
            "end1":        bins["end"].iloc[i],
            "start2":      bins["start"].iloc[j],
            "end2":        bins["end"].iloc[j],
            "enrichment":  round(float(enr), 3),
            "pvalue":      float(pv),
        })

    return pd.DataFrame(loops)

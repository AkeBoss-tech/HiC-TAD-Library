"""Chromatin loop detection using a donut background model.

Implements the core idea of Rao et al. (2014) / HICCUPS:
    1. For each pixel (i, j) in the upper triangle at a valid genomic
       distance, compute a 'donut' local background — the mean of pixels
       in an annular neighbourhood that avoids the short-range diagonal
       (which reflects TAD-internal enrichment, not loops).
    2. Compute enrichment = observed / donut_background.
    3. Assign a z-score using the local background distribution.
    4. Apply Benjamini–Hochberg FDR correction across all candidate pixels.
    5. Merge adjacent enriched pixels into single loop calls.
"""

import cooler
import numpy as np
import pandas as pd
from scipy import ndimage, stats
from typing import Optional

from .tad_boundaries import parse_coordinates


# ---------------------------------------------------------------------------
# Background model
# ---------------------------------------------------------------------------

def _make_donut_kernel(inner_radius: int, outer_radius: int) -> np.ndarray:
    """Return a normalised 2-D donut kernel (ring-shaped average filter)."""
    size = 2 * outer_radius + 1
    kernel = np.zeros((size, size), dtype=float)
    center = outer_radius
    for i in range(size):
        for j in range(size):
            dist = max(abs(i - center), abs(j - center))
            if inner_radius < dist <= outer_radius:
                kernel[i, j] = 1.0
    total = kernel.sum()
    if total == 0:
        raise ValueError("Donut kernel has zero area; check inner/outer radius.")
    return kernel / total


def _donut_background(
    matrix: np.ndarray,
    inner_radius: int = 3,
    outer_radius: int = 7,
) -> np.ndarray:
    """Estimate per-pixel local background via a donut convolution."""
    kernel = _make_donut_kernel(inner_radius, outer_radius)
    mat = np.nan_to_num(matrix, nan=0.0)
    return ndimage.convolve(mat, kernel, mode='reflect')


# ---------------------------------------------------------------------------
# Main loop caller
# ---------------------------------------------------------------------------

def call_loops(
    clr: cooler.Cooler,
    region: str,
    min_dist_bp: int = 100_000,
    max_dist_bp: int = 2_000_000,
    enrichment_threshold: float = 1.5,
    zscore_threshold: float = 2.0,
    fdr_threshold: float = 0.1,
    inner_radius: int = 3,
    outer_radius: int = 7,
    merge_distance_bins: int = 2,
) -> pd.DataFrame:
    """Call chromatin loops in a genomic region.

    Parameters
    ----------
    clr : cooler.Cooler
    region : str
    min_dist_bp : int     minimum loop span (bp)
    max_dist_bp : int     maximum loop span (bp)
    enrichment_threshold : float  O/E enrichment over donut background
    zscore_threshold : float      z-score threshold for significance
    fdr_threshold : float         maximum BH-adjusted p-value
    inner_radius : int            donut inner radius (bins)
    outer_radius : int            donut outer radius (bins)
    merge_distance_bins : int     merge peaks within this bin distance

    Returns
    -------
    pd.DataFrame
        chrom1, start1, end1, chrom2, start2, end2,
        observed, expected_bg, enrichment, zscore, pvalue, qvalue
    """
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    min_dist_bins = max(1, min_dist_bp // resolution)
    max_dist_bins = max_dist_bp // resolution

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)
    n = matrix.shape[0]

    bg = _donut_background(matrix, inner_radius=inner_radius, outer_radius=outer_radius)

    # ---- distance matrix (vectorised) ----
    ii, jj = np.indices((n, n))
    dist_mat = jj - ii   # signed; we want upper triangle so dist > 0

    valid_mask = (
        (dist_mat >= min_dist_bins) &
        (dist_mat <= max_dist_bins) &
        (bg > 0) &
        (matrix > 0)
    )

    obs_vals = matrix[valid_mask]
    bg_vals = bg[valid_mask]
    enr_vals = obs_vals / bg_vals

    # z-score against local background distribution
    bg_std = np.std(bg_vals)
    bg_mean = np.mean(bg_vals)
    if bg_std > 0:
        z_vals = (obs_vals - bg_mean) / bg_std
    else:
        z_vals = np.zeros_like(obs_vals)

    # Two-tailed normal p-values on z-scores
    pvals = stats.norm.sf(z_vals)

    # ---- BH FDR correction ----
    m = len(pvals)
    order = np.argsort(pvals)
    qvals = np.ones(m)
    cummin = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        qvals[idx] = min(cummin, pvals[idx] * m / (rank + 1))
        cummin = qvals[idx]

    # ---- Filter ----
    keep = (enr_vals > enrichment_threshold) & (z_vals > zscore_threshold) & (qvals <= fdr_threshold)
    if not keep.any():
        return pd.DataFrame(columns=[
            'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
            'observed', 'expected_bg', 'enrichment', 'zscore', 'pvalue', 'qvalue',
        ])

    ci_all, cj_all = np.where(valid_mask)
    ci = ci_all[keep]
    cj = cj_all[keep]
    enr_keep = enr_vals[keep]
    z_keep = z_vals[keep]
    p_keep = pvals[keep]
    q_keep = qvals[keep]

    # ---- Greedy merge of nearby pixels ----
    sort_idx = np.argsort(-enr_keep)
    used = np.zeros(len(ci), dtype=bool)
    loops = []

    for idx in sort_idx:
        if used[idx]:
            continue
        i, j = ci[idx], cj[idx]
        nearby = (
            (np.abs(ci - i) <= merge_distance_bins) &
            (np.abs(cj - j) <= merge_distance_bins)
        )
        used[nearby] = True
        loops.append({
            'chrom1': chrom,
            'start1': int(start + i * resolution),
            'end1': int(start + (i + 1) * resolution),
            'chrom2': chrom,
            'start2': int(start + j * resolution),
            'end2': int(start + (j + 1) * resolution),
            'observed': float(matrix[i, j]),
            'expected_bg': float(bg[i, j]),
            'enrichment': float(enr_keep[idx]),
            'zscore': float(z_keep[idx]),
            'pvalue': float(p_keep[idx]),
            'qvalue': float(q_keep[idx]),
        })

    return pd.DataFrame(loops)


# ---------------------------------------------------------------------------
# Pileup (aggregate peak analysis)
# ---------------------------------------------------------------------------

def loop_pileup(
    clr: cooler.Cooler,
    loops_df: pd.DataFrame,
    flank_bins: int = 10,
) -> np.ndarray:
    """Aggregate contact sub-matrices centred on loop anchors.

    Parameters
    ----------
    clr : cooler.Cooler
    loops_df : pd.DataFrame  output of call_loops
    flank_bins : int  bins to show on each side of each anchor

    Returns
    -------
    ndarray  shape (2*flank_bins+1, 2*flank_bins+1) — averaged pileup
    """
    resolution = clr.binsize
    size = 2 * flank_bins + 1
    stack = []

    for _, row in loops_df.iterrows():
        chrom = row['chrom1']
        ci = (row['start1'] + resolution // 2) // resolution
        cj = (row['start2'] + resolution // 2) // resolution
        chrom_bins = clr.chromsizes[chrom] // resolution

        lo_i, hi_i = ci - flank_bins, ci + flank_bins + 1
        lo_j, hi_j = cj - flank_bins, cj + flank_bins + 1

        if lo_i < 0 or hi_i > chrom_bins or lo_j < 0 or hi_j > chrom_bins:
            continue

        snippet = clr.matrix(balance=True).fetch(
            f"{chrom}:{lo_i * resolution}-{hi_i * resolution}",
            f"{chrom}:{lo_j * resolution}-{hi_j * resolution}",
        )
        if snippet.shape == (size, size):
            stack.append(np.nan_to_num(snippet, nan=0.0))

    if not stack:
        return np.zeros((size, size))
    return np.mean(stack, axis=0)

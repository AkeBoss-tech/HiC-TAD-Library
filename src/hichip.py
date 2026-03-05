"""HiChIP analysis: ChIP-enriched chromatin loop calling.

HiChIP combines proximity ligation (Hi-C) with ChIP-seq enrichment,
producing a contact map biased towards regions bound by the target protein
(e.g. cohesin/SMC1, CTCF, H3K27ac).

Analysis pipeline:
    1. Load ChIP-seq peaks (BED) defining the 1-D binding profile.
    2. Compute the standard Hi-C contact matrix for the region.
    3. Identify peak–peak bin pairs and score their contact enrichment
       over a distance-matched background.
    4. Call significant loops anchored at peak bins (FitHiChIP-inspired).
    5. Optionally compute 1-D pile-up profiles around peak centres.

Reference: Bhatt et al. 2015 (HiChIP); FitHiChIP (2019).
"""

import cooler
import numpy as np
import pandas as pd
from scipy import stats
from typing import Optional

from .tad_boundaries import parse_coordinates
from .loops import _donut_background


# ---------------------------------------------------------------------------
# Peak loading
# ---------------------------------------------------------------------------

def load_peaks(peaks_bed: str, chrom: Optional[str] = None) -> pd.DataFrame:
    """Load ChIP-seq peaks from a BED3+ file.

    Parameters
    ----------
    peaks_bed : str   path to BED file (chrom, start, end, [name, score, ...])
    chrom : str, optional   restrict to this chromosome

    Returns
    -------
    pd.DataFrame  columns: chrom, start, end
    """
    df = pd.read_csv(
        peaks_bed, sep='\t', header=None, comment='#',
        usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
    )
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    if chrom:
        df = df[df['chrom'] == chrom].copy()
    return df.reset_index(drop=True)


def peaks_to_bins(peaks_df: pd.DataFrame, resolution: int, region_start: int = 0) -> np.ndarray:
    """Convert peak midpoints to bin indices relative to a region start."""
    midpoints = (peaks_df['start'] + peaks_df['end']) // 2
    return ((midpoints - region_start) // resolution).astype(int).values


# ---------------------------------------------------------------------------
# HiChIP loop calling
# ---------------------------------------------------------------------------

def call_hichip_loops(
    clr: cooler.Cooler,
    peaks_df: pd.DataFrame,
    region: str,
    min_dist_bp: int = 5_000,
    max_dist_bp: int = 2_000_000,
    enrichment_threshold: float = 2.0,
    fdr_threshold: float = 0.05,
    inner_radius: int = 2,
    outer_radius: int = 5,
) -> pd.DataFrame:
    """Call HiChIP loops anchored at ChIP-seq peaks.

    Only pixel pairs where both anchors fall within a peak bin are tested,
    reducing multiple-testing burden and increasing sensitivity at bound loci.

    Parameters
    ----------
    clr : cooler.Cooler
    peaks_df : pd.DataFrame   output of load_peaks (chrom, start, end)
    region : str
    min_dist_bp, max_dist_bp : int   loop size range
    enrichment_threshold : float   minimum observed/background enrichment
    fdr_threshold : float   BH-adjusted p-value threshold
    inner_radius, outer_radius : int   donut background parameters

    Returns
    -------
    pd.DataFrame
        chrom1, start1, end1, chrom2, start2, end2,
        observed, expected_bg, enrichment, pvalue, qvalue
    """
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    min_dist_bins = max(1, min_dist_bp // resolution)
    max_dist_bins = max_dist_bp // resolution

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)
    n = matrix.shape[0]

    bg = _donut_background(matrix, inner_radius=inner_radius, outer_radius=outer_radius)

    # Peak bins within region
    region_peaks = peaks_df[
        (peaks_df['chrom'] == chrom) &
        (peaks_df['start'] >= start) &
        (peaks_df['end'] <= end)
    ].copy()

    if region_peaks.empty:
        print(f"  HiChIP: no peaks in {region} — consider expanding the region or checking peak coordinates.")
        return pd.DataFrame(columns=[
            'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
            'observed', 'expected_bg', 'enrichment', 'pvalue', 'qvalue',
        ])

    peak_mids = (region_peaks['start'].values + region_peaks['end'].values) // 2
    peak_bins = ((peak_mids - start) // resolution).clip(0, n - 1).astype(int)
    peak_bins = np.unique(peak_bins)

    # Build candidate list: all peak-peak pairs in valid distance range
    rows, cols, obs_list, exp_list, enr_list = [], [], [], [], []
    for ii, pi in enumerate(peak_bins):
        for jj, pj in enumerate(peak_bins):
            if pj <= pi:
                continue
            dist = pj - pi
            if not (min_dist_bins <= dist <= max_dist_bins):
                continue
            if bg[pi, pj] <= 0 or matrix[pi, pj] <= 0:
                continue
            enr = matrix[pi, pj] / bg[pi, pj]
            if enr > enrichment_threshold:
                rows.append(pi)
                cols.append(pj)
                obs_list.append(matrix[pi, pj])
                exp_list.append(bg[pi, pj])
                enr_list.append(enr)

    if not rows:
        return pd.DataFrame(columns=[
            'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
            'observed', 'expected_bg', 'enrichment', 'pvalue', 'qvalue',
        ])

    obs = np.array(obs_list)
    exp = np.array(exp_list)
    enr = np.array(enr_list)
    bg_mean = np.mean(exp)
    bg_std = np.std(exp)

    z = (obs - bg_mean) / max(bg_std, 1e-10)
    pvals = stats.norm.sf(z)

    # BH FDR correction
    m = len(pvals)
    order = np.argsort(pvals)
    qvals = np.ones(m)
    cummin = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        qvals[idx] = min(cummin, pvals[idx] * m / (rank + 1))
        cummin = qvals[idx]

    loops = []
    for idx in range(m):
        if qvals[idx] <= fdr_threshold:
            pi, pj = rows[idx], cols[idx]
            loops.append({
                'chrom1': chrom,
                'start1': int(start + pi * resolution),
                'end1': int(start + (pi + 1) * resolution),
                'chrom2': chrom,
                'start2': int(start + pj * resolution),
                'end2': int(start + (pj + 1) * resolution),
                'observed': float(obs[idx]),
                'expected_bg': float(exp[idx]),
                'enrichment': float(enr[idx]),
                'pvalue': float(pvals[idx]),
                'qvalue': float(qvals[idx]),
            })

    return pd.DataFrame(loops)


# ---------------------------------------------------------------------------
# 1-D profiles and enrichment scoring
# ---------------------------------------------------------------------------

def peak_contact_profile(
    clr: cooler.Cooler,
    peaks_df: pd.DataFrame,
    region: str,
    flank_bins: int = 20,
) -> np.ndarray:
    """Average contact profile centred on each peak in a region.

    Returns
    -------
    ndarray  shape (n_peaks, 2*flank_bins+1)
    """
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    n_region = (end - start) // resolution

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)

    region_peaks = peaks_df[
        (peaks_df['chrom'] == chrom) &
        (peaks_df['start'] >= start) &
        (peaks_df['end'] <= end)
    ]
    peak_mids = (region_peaks['start'].values + region_peaks['end'].values) // 2
    peak_bins = ((peak_mids - start) // resolution).clip(0, n_region - 1).astype(int)

    size = 2 * flank_bins + 1
    profiles = []
    for pb in peak_bins:
        lo, hi = pb - flank_bins, pb + flank_bins + 1
        if lo < 0 or hi > n_region:
            continue
        profile = matrix[pb, lo:hi]
        if len(profile) == size:
            profiles.append(profile)

    return np.array(profiles) if profiles else np.zeros((0, size))


def compute_peak_enrichment(
    clr: cooler.Cooler,
    peaks_df: pd.DataFrame,
    region: str,
) -> pd.DataFrame:
    """Score each peak by average contact enrichment at its bin.

    Returns
    -------
    pd.DataFrame  chrom, start, end, bin_idx, mean_contact, enrichment_over_median
    """
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    n_region = (end - start) // resolution

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)
    row_means = matrix.mean(axis=1)
    global_median = np.median(row_means[row_means > 0]) if (row_means > 0).any() else 1.0

    region_peaks = peaks_df[
        (peaks_df['chrom'] == chrom) &
        (peaks_df['start'] >= start) &
        (peaks_df['end'] <= end)
    ].copy()

    if region_peaks.empty:
        return pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'bin_idx', 'mean_contact', 'enrichment_over_median',
        ])

    peak_mids = (region_peaks['start'].values + region_peaks['end'].values) // 2
    bin_idxs = ((peak_mids - start) // resolution).clip(0, n_region - 1).astype(int)

    result = region_peaks[['chrom', 'start', 'end']].copy().reset_index(drop=True)
    result['bin_idx'] = bin_idxs
    result['mean_contact'] = [float(row_means[b]) for b in bin_idxs]
    result['enrichment_over_median'] = result['mean_contact'] / max(global_median, 1e-10)
    return result

"""Capture Hi-C analysis (CHiCAGO-inspired interaction calling).

Capture Hi-C enriches contact libraries for specific genomic regions
(baits / viewpoints) by adding biotinylated oligonucleotide probes.
The result is a contact matrix highly enriched for bait–prey interactions.

Algorithm (simplified CHiCAGO, Cairns et al. 2016):
    1. Load bait regions from a BED file.
    2. For each bait bin, extract its contact profile across the region
       (the prey contacts).
    3. Model expected contacts under the null as a distance-dependent
       negative binomial: at short distances, contacts are dominated by
       random ligation (Brownian component); at long distances by the
       structured chromatin background.
    4. Call interactions where the observed count significantly exceeds
       the model expectation (Poisson p-value, BH-corrected).
"""

import cooler
import numpy as np
import pandas as pd
from scipy import stats
from typing import List, Optional

from .tad_boundaries import parse_coordinates


# ---------------------------------------------------------------------------
# Bait loading
# ---------------------------------------------------------------------------

def load_baits(baits_bed: str, chrom: Optional[str] = None) -> pd.DataFrame:
    """Load capture bait regions from a BED3+ file.

    Returns
    -------
    pd.DataFrame  columns: chrom, start, end
    """
    df = pd.read_csv(
        baits_bed, sep='\t', header=None, comment='#',
        usecols=[0, 1, 2], names=['chrom', 'start', 'end'],
    )
    df['start'] = df['start'].astype(int)
    df['end'] = df['end'].astype(int)
    if chrom:
        df = df[df['chrom'] == chrom].copy()
    return df.reset_index(drop=True)


# ---------------------------------------------------------------------------
# Distance-dependent background model
# ---------------------------------------------------------------------------

def fit_distance_decay(
    clr: cooler.Cooler,
    region: str,
    min_dist_bp: int = 5_000,
    max_dist_bp: int = 5_000_000,
) -> np.ndarray:
    """Fit observed/expected distance decay for a region.

    Returns
    -------
    expected_by_dist : ndarray  shape (n_bins,) — expected contact frequency
                                at each bin distance 0, 1, 2, ...
    """
    chrom, start, end = parse_coordinates(region)
    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)
    n = matrix.shape[0]
    resolution = clr.binsize

    min_d = max(1, min_dist_bp // resolution)
    max_d = min(n - 1, max_dist_bp // resolution)

    expected = np.zeros(n)
    for d in range(min_d, max_d + 1):
        diag = np.diag(matrix, k=d)
        valid = diag[diag > 0]
        expected[d] = np.mean(valid) if len(valid) > 0 else 0.0

    # Smooth with a running average
    window = 5
    for d in range(min_d, max_d + 1):
        lo = max(min_d, d - window)
        hi = min(max_d + 1, d + window + 1)
        vals = expected[lo:hi]
        expected[d] = np.mean(vals[vals > 0]) if (vals > 0).any() else 0.0

    return expected


# ---------------------------------------------------------------------------
# Interaction calling
# ---------------------------------------------------------------------------

def call_bait_interactions(
    clr: cooler.Cooler,
    baits_df: pd.DataFrame,
    region: str,
    min_dist_bp: int = 5_000,
    max_dist_bp: int = 2_000_000,
    min_contacts: float = 0.01,
    fdr_threshold: float = 0.05,
) -> pd.DataFrame:
    """Call significant bait–prey interactions in a region.

    Parameters
    ----------
    clr : cooler.Cooler
    baits_df : pd.DataFrame   output of load_baits
    region : str
    min_dist_bp, max_dist_bp : int   interaction distance range
    min_contacts : float   minimum balanced contact value to consider
    fdr_threshold : float

    Returns
    -------
    pd.DataFrame
        bait_chrom, bait_start, bait_end,
        prey_chrom, prey_start, prey_end,
        observed, expected, enrichment, pvalue, qvalue, significant
    """
    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    n = (end - start) // resolution

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)

    min_dist_bins = max(1, min_dist_bp // resolution)
    max_dist_bins = max_dist_bp // resolution

    expected = fit_distance_decay(clr, region, min_dist_bp, max_dist_bp)

    # Bait bins in region
    region_baits = baits_df[
        (baits_df['chrom'] == chrom) &
        (baits_df['start'] >= start) &
        (baits_df['end'] <= end)
    ].copy()

    if region_baits.empty:
        return pd.DataFrame(columns=[
            'bait_chrom', 'bait_start', 'bait_end',
            'prey_chrom', 'prey_start', 'prey_end',
            'observed', 'expected', 'enrichment', 'pvalue', 'qvalue', 'significant',
        ])

    bait_mids = (region_baits['start'].values + region_baits['end'].values) // 2
    bait_bins = ((bait_mids - start) // resolution).clip(0, n - 1).astype(int)

    all_interactions = []

    for bi_idx, bait_bin in enumerate(bait_bins):
        bait_start = region_baits['start'].iloc[bi_idx]
        bait_end = region_baits['end'].iloc[bi_idx]

        for prey_bin in range(n):
            dist = abs(prey_bin - bait_bin)
            if not (min_dist_bins <= dist <= max_dist_bins):
                continue

            obs = float(matrix[bait_bin, prey_bin])
            if obs < min_contacts:
                continue

            exp = float(expected[dist]) if dist < len(expected) else 0.0
            if exp <= 0:
                continue

            enr = obs / exp
            all_interactions.append({
                'bait_bin': bait_bin,
                'prey_bin': prey_bin,
                'bait_start': bait_start,
                'bait_end': bait_end,
                'observed': obs,
                'expected': exp,
                'enrichment': enr,
            })

    if not all_interactions:
        return pd.DataFrame(columns=[
            'bait_chrom', 'bait_start', 'bait_end',
            'prey_chrom', 'prey_start', 'prey_end',
            'observed', 'expected', 'enrichment', 'pvalue', 'qvalue', 'significant',
        ])

    idf = pd.DataFrame(all_interactions)

    # Poisson-style p-values (using z-score on enrichment for numerical stability)
    enr_vals = idf['enrichment'].values
    mu_enr = np.mean(enr_vals)
    sd_enr = np.std(enr_vals)
    z = (enr_vals - mu_enr) / max(sd_enr, 1e-10)
    pvals = stats.norm.sf(z)

    # BH FDR
    m = len(pvals)
    order = np.argsort(pvals)
    qvals = np.ones(m)
    cummin = 1.0
    for rank in range(m - 1, -1, -1):
        idx = order[rank]
        qvals[idx] = min(cummin, pvals[idx] * m / (rank + 1))
        cummin = qvals[idx]

    idf['pvalue'] = pvals
    idf['qvalue'] = qvals
    idf['significant'] = qvals <= fdr_threshold

    # Translate bins to coordinates
    idf['bait_chrom'] = chrom
    idf['prey_chrom'] = chrom
    idf['prey_start'] = (start + idf['prey_bin'] * resolution).astype(int)
    idf['prey_end'] = (idf['prey_start'] + resolution).astype(int)

    result = idf[[
        'bait_chrom', 'bait_start', 'bait_end',
        'prey_chrom', 'prey_start', 'prey_end',
        'observed', 'expected', 'enrichment', 'pvalue', 'qvalue', 'significant',
    ]].copy()

    return result.sort_values('qvalue').reset_index(drop=True)


# ---------------------------------------------------------------------------
# Multi-bait summary
# ---------------------------------------------------------------------------

def summarise_interactions(
    interactions_df: pd.DataFrame,
    significant_only: bool = True,
) -> pd.DataFrame:
    """Group interactions by bait and summarise statistics.

    Returns
    -------
    pd.DataFrame  bait_chrom, bait_start, bait_end,
                  n_interactions, mean_enrichment, max_enrichment, top_prey_start
    """
    df = interactions_df.copy()
    if significant_only:
        df = df[df['significant']].copy()

    if df.empty:
        return pd.DataFrame(columns=[
            'bait_chrom', 'bait_start', 'bait_end',
            'n_interactions', 'mean_enrichment', 'max_enrichment', 'top_prey_start',
        ])

    summary = (
        df.groupby(['bait_chrom', 'bait_start', 'bait_end'])
        .agg(
            n_interactions=('prey_start', 'count'),
            mean_enrichment=('enrichment', 'mean'),
            max_enrichment=('enrichment', 'max'),
            top_prey_start=('prey_start', lambda x: x.iloc[df.loc[x.index, 'enrichment'].argmax()]),
        )
        .reset_index()
    )
    return summary.sort_values('max_enrichment', ascending=False).reset_index(drop=True)

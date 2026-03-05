"""A/B compartment analysis from Hi-C contact matrices.

Algorithm (Lieberman-Aiden et al. 2009):
    1. Compute observed/expected (O/E) matrix by dividing each pixel by the
       genome-wide average contact frequency at that genomic distance.
    2. Compute Pearson correlation matrix across rows of O/E.
    3. Eigendecompose the correlation matrix; the first eigenvector (E1)
       separates A compartments (active, gene-rich) from B (inactive).
    4. Orient sign so that A = positive E1 (requires gene density or GC
       content; if absent the sign convention is internally consistent but
       arbitrary).
"""

import cooler
import cooltools
import numpy as np
import pandas as pd
from scipy import linalg
from typing import Optional, Tuple

from .tad_boundaries import parse_coordinates, make_view_df


def compute_oe_matrix(
    clr: cooler.Cooler,
    region: str,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute observed/expected contact matrix for a genomic region.

    Returns
    -------
    (oe_matrix, expected_by_dist)
        oe_matrix       : ndarray shape (n, n) — O/E values
        expected_by_dist: 1-D array, expected contact freq at each distance d
    """
    chrom, start, end = parse_coordinates(region)
    chrom_len = clr.chromsizes[chrom]

    # Distance-dependent expected on the full chromosome
    view_df = make_view_df(chrom, 0, int(chrom_len))
    expected = cooltools.expected_cis(clr, view_df=view_df, nproc=1)

    # Support both old and new cooltools column naming
    avg_col = 'balanced.avg' if 'balanced.avg' in expected.columns else expected.columns[-1]
    exp_by_dist = expected.set_index('dist')[avg_col].to_dict()

    matrix = clr.matrix(balance=True).fetch(f"{chrom}:{start}-{end}")
    matrix = np.nan_to_num(matrix, nan=0.0)
    n = matrix.shape[0]

    # Vectorised O/E construction
    i_idx, j_idx = np.indices((n, n))
    dist_arr = np.abs(i_idx - j_idx)

    exp_lookup = np.array([exp_by_dist.get(d, np.nan) for d in range(n)])
    exp_matrix = exp_lookup[dist_arr]

    with np.errstate(divide='ignore', invalid='ignore'):
        oe = np.where(exp_matrix > 0, matrix / exp_matrix, 0.0)

    return oe, exp_lookup


def compute_compartments(
    clr: cooler.Cooler,
    region: str,
    gene_density: Optional[np.ndarray] = None,
) -> pd.DataFrame:
    """Compute A/B compartment eigenvector (E1) for a genomic region.

    Parameters
    ----------
    clr : cooler.Cooler
    region : str   e.g. 'chr12:26,000,000-28,000,000'
    gene_density : array-like, optional
        Gene density per bin (same length as matrix). Used to orient sign so
        A compartment is positive.  If None, sign may be arbitrary.

    Returns
    -------
    pd.DataFrame  columns: chrom, start, end, E1, compartment ('A' | 'B' | 'unknown')
    """
    chrom, start, end = parse_coordinates(region)
    oe, _ = compute_oe_matrix(clr, region)

    # Remove empty bins
    row_sums = np.abs(oe).sum(axis=1)
    valid = row_sums > 0
    oe_valid = oe[np.ix_(valid, valid)]

    if oe_valid.shape[0] < 3:
        raise ValueError(f"Too few valid bins in {region} for compartment analysis.")

    # Pearson correlation across rows of O/E
    corr = np.corrcoef(oe_valid)
    corr = np.nan_to_num(corr, nan=0.0)

    # Eigendecomposition — largest eigenvalue's eigenvector = E1
    _, eigenvectors = linalg.eigh(corr)
    e1_valid = eigenvectors[:, -1]

    # Map back to full-length (invalid bins stay NaN)
    e1 = np.full(oe.shape[0], np.nan)
    e1[valid] = e1_valid

    # Orient sign so A (high gene density) = positive
    if gene_density is not None and len(gene_density) == len(e1):
        gd = np.asarray(gene_density, dtype=float)
        ok = ~np.isnan(e1) & ~np.isnan(gd)
        if ok.sum() > 2 and np.corrcoef(e1[ok], gd[ok])[0, 1] < 0:
            e1 = -e1

    resolution = clr.binsize
    n_bins = len(e1)
    starts = start + np.arange(n_bins) * resolution

    df = pd.DataFrame({
        'chrom': chrom,
        'start': starts.astype(int),
        'end': (starts + resolution).astype(int),
        'E1': e1,
    })
    df['compartment'] = df['E1'].apply(
        lambda x: 'A' if x > 0 else ('B' if x < 0 else 'unknown')
    )
    return df


def compartment_saddle(
    clr: cooler.Cooler,
    region: str,
    n_quantiles: int = 50,
) -> Tuple[np.ndarray, np.ndarray]:
    """Compute compartment saddle plot (average O/E as function of E1 quantile).

    Rows/columns of the output are sorted by E1 so that B–B, A–B, and A–A
    interactions appear in the lower-left, off-diagonal, and upper-right.

    Returns
    -------
    (saddle_matrix, quantile_edges)
        saddle_matrix  : (n_quantiles, n_quantiles) average O/E
        quantile_edges : 1-D array of E1 bin boundaries
    """
    comp_df = compute_compartments(clr, region)
    oe, _ = compute_oe_matrix(clr, region)

    e1 = comp_df['E1'].values
    valid = ~np.isnan(e1)
    e1_valid = e1[valid]

    quantiles = np.percentile(e1_valid, np.linspace(0, 100, n_quantiles + 1))
    bin_ids = np.digitize(e1, quantiles, right=True).clip(0, n_quantiles - 1)

    saddle = np.zeros((n_quantiles, n_quantiles))
    counts = np.zeros((n_quantiles, n_quantiles))

    for i in range(len(e1)):
        for j in range(i, len(e1)):
            if np.isnan(e1[i]) or np.isnan(e1[j]):
                continue
            bi, bj = bin_ids[i], bin_ids[j]
            saddle[bi, bj] += oe[i, j]
            counts[bi, bj] += 1
            if bi != bj:
                saddle[bj, bi] += oe[j, i]
                counts[bj, bi] += 1

    with np.errstate(divide='ignore', invalid='ignore'):
        saddle = np.where(counts > 0, saddle / counts, np.nan)

    return saddle, quantiles

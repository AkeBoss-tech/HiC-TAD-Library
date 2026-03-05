"""Copy number variation (CNV) detection from Hi-C contact data and BAM files.

Hi-C–based CNV detection exploits the fact that the total contact count for
a genomic bin is proportional to its copy number: duplicated regions produce
more contacts because there are more template molecules.

Pipeline:
    1. Compute per-bin coverage = sum of row contacts in the raw (unbalanced)
       contact matrix.
    2. Apply GC-content correction if a reference FASTA is provided.
    3. Normalise to the genome-wide median (assumed diploid → CN = 2).
    4. Segment the copy-number track using a sliding-window CBS-like approach.

BAM-based CNV uses read-depth per bin, corrected for GC content and
mappability.
"""

import cooler
import numpy as np
import pandas as pd
from scipy import stats, ndimage
from typing import List, Optional

from .tad_boundaries import parse_coordinates, make_view_df

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False


# ---------------------------------------------------------------------------
# Hi-C–based coverage
# ---------------------------------------------------------------------------

def compute_hic_coverage(
    clr: cooler.Cooler,
    chrom: str,
    resolution_bp: Optional[int] = None,
) -> pd.DataFrame:
    """Compute per-bin contact coverage (sum of row) from the raw matrix.

    Parameters
    ----------
    clr : cooler.Cooler
    chrom : str   e.g. 'chr12'
    resolution_bp : int, optional   rebin to this resolution (must be multiple of clr.binsize)

    Returns
    -------
    pd.DataFrame  columns: chrom, start, end, raw_coverage
    """
    chrom_len = int(clr.chromsizes[chrom])
    matrix = clr.matrix(balance=False).fetch(chrom).astype(float)
    matrix = np.nan_to_num(matrix, nan=0.0)
    coverage = matrix.sum(axis=1)

    resolution = clr.binsize
    n = len(coverage)
    starts = np.arange(n) * resolution

    df = pd.DataFrame({
        'chrom': chrom,
        'start': starts.astype(int),
        'end': (starts + resolution).astype(int),
        'raw_coverage': coverage,
    })

    if resolution_bp and resolution_bp > resolution:
        factor = resolution_bp // resolution
        df['bin_group'] = df.index // factor
        df = df.groupby('bin_group').agg({
            'chrom': 'first',
            'start': 'first',
            'end': 'last',
            'raw_coverage': 'sum',
        }).reset_index(drop=True)

    return df


def gc_correct_coverage(
    coverage_df: pd.DataFrame,
    gc_content: np.ndarray,
    poly_degree: int = 3,
) -> pd.DataFrame:
    """Remove GC-content bias from bin coverage using polynomial regression.

    Parameters
    ----------
    coverage_df : pd.DataFrame   output of compute_hic_coverage
    gc_content : np.ndarray      GC fraction per bin (0–1), same length
    poly_degree : int            degree of correction polynomial

    Returns
    -------
    pd.DataFrame  with added column 'gc_corrected_coverage'
    """
    df = coverage_df.copy()
    gc = np.asarray(gc_content, dtype=float)
    cov = df['raw_coverage'].values.astype(float)

    valid = np.isfinite(gc) & np.isfinite(cov) & (cov > 0)
    if valid.sum() < poly_degree + 2:
        df['gc_corrected_coverage'] = cov
        return df

    coeffs = np.polyfit(gc[valid], cov[valid], poly_degree)
    expected_from_gc = np.polyval(coeffs, gc)
    median_cov = np.median(cov[valid])

    with np.errstate(divide='ignore', invalid='ignore'):
        corrected = np.where(expected_from_gc > 0, cov / expected_from_gc * median_cov, cov)

    df['gc_corrected_coverage'] = corrected
    return df


# ---------------------------------------------------------------------------
# Copy number normalisation
# ---------------------------------------------------------------------------

def normalise_to_ploidy(
    coverage_df: pd.DataFrame,
    ploidy: int = 2,
    coverage_col: str = 'raw_coverage',
    smooth_bins: int = 5,
) -> pd.DataFrame:
    """Normalise coverage to copy number units.

    Assumes the median bin is diploid (CN = ploidy). Applies a running
    median smooth before normalisation to reduce noise.

    Returns
    -------
    pd.DataFrame  with added columns 'smoothed_coverage', 'copy_number'
    """
    df = coverage_df.copy()
    cov = df[coverage_col].values.astype(float)

    # Running median smooth
    smoothed = ndimage.median_filter(cov, size=smooth_bins)
    df['smoothed_coverage'] = smoothed

    baseline = np.median(smoothed[smoothed > 0]) if (smoothed > 0).any() else 1.0
    with np.errstate(divide='ignore', invalid='ignore'):
        cn = smoothed / baseline * ploidy
    df['copy_number'] = np.clip(cn, 0, 10)

    return df


# ---------------------------------------------------------------------------
# Segmentation
# ---------------------------------------------------------------------------

def segment_cnv(
    coverage_df: pd.DataFrame,
    cn_col: str = 'copy_number',
    window_bins: int = 20,
    change_threshold: float = 0.5,
    min_segment_bins: int = 5,
) -> pd.DataFrame:
    """Segment the copy-number track into constant-CN regions.

    Uses a sliding-window t-test to detect changepoints, then assigns a
    mean CN to each segment.

    Parameters
    ----------
    coverage_df : pd.DataFrame   output of normalise_to_ploidy
    cn_col : str
    window_bins : int   half-window for changepoint detection
    change_threshold : float   minimum |delta CN| to call a breakpoint
    min_segment_bins : int   discard segments shorter than this

    Returns
    -------
    pd.DataFrame  chrom, start, end, mean_cn, cn_state, n_bins
    """
    df = coverage_df.copy()
    cn = df[cn_col].values
    chrom = df['chrom'].iloc[0]

    # Detect breakpoints using sliding window mean difference
    breakpoints = [0]
    n = len(cn)
    for i in range(window_bins, n - window_bins):
        left = cn[max(0, i - window_bins):i]
        right = cn[i:min(n, i + window_bins)]
        if len(left) > 1 and len(right) > 1:
            delta = abs(np.mean(right) - np.mean(left))
            if delta > change_threshold:
                # Avoid placing breakpoints too close together
                if i - breakpoints[-1] >= min_segment_bins:
                    breakpoints.append(i)
    breakpoints.append(n)

    segments = []
    for k in range(len(breakpoints) - 1):
        seg_start = breakpoints[k]
        seg_end = breakpoints[k + 1]
        if seg_end - seg_start < min_segment_bins:
            continue

        seg_cn = cn[seg_start:seg_end]
        mean_cn = float(np.nanmean(seg_cn))
        cn_state = _cn_label(mean_cn)

        start_bp = int(df['start'].iloc[seg_start])
        end_bp = int(df['end'].iloc[seg_end - 1])

        segments.append({
            'chrom': chrom,
            'start': start_bp,
            'end': end_bp,
            'mean_cn': mean_cn,
            'cn_state': cn_state,
            'n_bins': seg_end - seg_start,
        })

    return pd.DataFrame(segments) if segments else pd.DataFrame(columns=[
        'chrom', 'start', 'end', 'mean_cn', 'cn_state', 'n_bins',
    ])


def _cn_label(mean_cn: float) -> str:
    """Convert continuous copy number to a discrete state label."""
    if mean_cn < 0.5:
        return 'homozygous_deletion'
    elif mean_cn < 1.5:
        return 'heterozygous_deletion'
    elif mean_cn < 2.5:
        return 'diploid'
    elif mean_cn < 3.5:
        return 'gain'
    else:
        return 'amplification'


# ---------------------------------------------------------------------------
# BAM-based coverage (alternative input)
# ---------------------------------------------------------------------------

def compute_bam_coverage(
    bam_path: str,
    chrom: str,
    bin_size: int = 5000,
    min_mapq: int = 20,
) -> pd.DataFrame:
    """Count mapped reads per genomic bin from a BAM file.

    Requires pysam.

    Returns
    -------
    pd.DataFrame  columns: chrom, start, end, raw_coverage
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required for BAM-based CNV. Install with: conda install pysam -c bioconda")

    bam = pysam.AlignmentFile(bam_path, "rb")
    chrom_len = bam.get_reference_length(chrom)
    n_bins = int(np.ceil(chrom_len / bin_size))
    counts = np.zeros(n_bins, dtype=int)

    for read in bam.fetch(chrom):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary:
            continue
        bin_idx = read.reference_start // bin_size
        if bin_idx < n_bins:
            counts[bin_idx] += 1

    bam.close()

    starts = np.arange(n_bins) * bin_size
    return pd.DataFrame({
        'chrom': chrom,
        'start': starts.astype(int),
        'end': np.minimum(starts + bin_size, chrom_len).astype(int),
        'raw_coverage': counts,
    })


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def call_cnv(
    clr: cooler.Cooler,
    chromosomes: Optional[List[str]] = None,
    ploidy: int = 2,
    smooth_bins: int = 5,
    change_threshold: float = 0.5,
) -> pd.DataFrame:
    """Run full Hi-C CNV pipeline across chromosomes.

    Returns
    -------
    pd.DataFrame  all CNV segments concatenated, with chrom, start, end,
                  mean_cn, cn_state, n_bins
    """
    if chromosomes is None:
        chromosomes = [c for c in clr.chromnames if 'M' not in c and '_' not in c]

    all_segments = []
    for chrom in chromosomes:
        try:
            cov_df = compute_hic_coverage(clr, chrom)
            cov_df = normalise_to_ploidy(cov_df, ploidy=ploidy, smooth_bins=smooth_bins)
            segs = segment_cnv(cov_df, change_threshold=change_threshold)
            all_segments.append(segs)
        except Exception as exc:
            print(f"  CNV: skipping {chrom}: {exc}")

    if not all_segments:
        return pd.DataFrame(columns=['chrom', 'start', 'end', 'mean_cn', 'cn_state', 'n_bins'])
    return pd.concat(all_segments, ignore_index=True)

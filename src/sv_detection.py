"""Structural variant detection from Hi-C contact data and BAM files.

Two detection strategies are implemented:

Hi-C–based (works from .mcool):
    - Translocations: inter-chromosomal contacts between two regions that
      exceed a z-score threshold relative to the chromosome-pair background.
    - Large deletions: contiguous runs of near-zero coverage bins flanked by
      elevated contacts bridging the gap.
    - Inversions: directional index sign-flip signature within a TAD.

BAM–based (requires pysam):
    - Discordant read pairs (unexpected orientation or insert size).
    - Split reads (large soft-clipped alignments spanning breakpoints).
"""

import cooler
import numpy as np
import pandas as pd
from scipy import stats
from typing import List, Optional, Tuple

from .tad_boundaries import parse_coordinates, make_view_df

# pysam is optional — only needed for BAM-based detection
try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False


# ---------------------------------------------------------------------------
# Hi-C based: translocation detection
# ---------------------------------------------------------------------------

def detect_translocations_hic(
    clr: cooler.Cooler,
    chromosomes: Optional[List[str]] = None,
    resolution_bp: Optional[int] = None,
    zscore_threshold: float = 5.0,
    min_contacts: int = 3,
) -> pd.DataFrame:
    """Identify putative translocations from inter-chromosomal Hi-C contacts.

    High-confidence translocations produce a focal cluster of inter-
    chromosomal contacts between two genomic loci. We z-score each
    chromosome-pair contact matrix and flag (chrom1, bin1, chrom2, bin2)
    combinations that exceed the threshold.

    Parameters
    ----------
    clr : cooler.Cooler
    chromosomes : list of str, optional   chromosomes to scan (default: all)
    resolution_bp : int, optional   coarser resolution to fetch (default: clr.binsize)
    zscore_threshold : float
    min_contacts : int   minimum raw contact count to consider

    Returns
    -------
    pd.DataFrame
        chrom1, start1, end1, chrom2, start2, end2, contacts, zscore, sv_type
    """
    if chromosomes is None:
        chromosomes = [c for c in clr.chromnames if 'M' not in c and '_' not in c]

    resolution = resolution_bp or clr.binsize
    results = []

    for i, c1 in enumerate(chromosomes):
        for j, c2 in enumerate(chromosomes):
            if j <= i:
                continue
            try:
                mat = clr.matrix(balance=False).fetch(c1, c2).astype(float)
            except Exception:
                continue

            if mat.size == 0:
                continue

            # Rebin if requested
            if resolution != clr.binsize and resolution > clr.binsize:
                factor = resolution // clr.binsize
                n1 = (mat.shape[0] // factor) * factor
                n2 = (mat.shape[1] // factor) * factor
                mat = mat[:n1, :n2]
                mat = mat.reshape(n1 // factor, factor, n2 // factor, factor).sum(axis=(1, 3))

            flat = mat.flatten()
            mu = flat.mean()
            sigma = flat.std()
            if sigma == 0:
                continue

            z_mat = (mat - mu) / sigma
            hot = np.argwhere(z_mat > zscore_threshold)

            for bi, bj in hot:
                if mat[bi, bj] >= min_contacts:
                    results.append({
                        'chrom1': c1,
                        'start1': int(bi * resolution),
                        'end1': int((bi + 1) * resolution),
                        'chrom2': c2,
                        'start2': int(bj * resolution),
                        'end2': int((bj + 1) * resolution),
                        'contacts': float(mat[bi, bj]),
                        'zscore': float(z_mat[bi, bj]),
                        'sv_type': 'translocation',
                    })

    return pd.DataFrame(results) if results else pd.DataFrame(columns=[
        'chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2',
        'contacts', 'zscore', 'sv_type',
    ])


# ---------------------------------------------------------------------------
# Hi-C based: large deletion detection
# ---------------------------------------------------------------------------

def detect_deletions_hic(
    clr: cooler.Cooler,
    region: str,
    zero_coverage_threshold: float = 0.05,
    min_gap_bins: int = 2,
    flank_enrichment_threshold: float = 1.5,
) -> pd.DataFrame:
    """Identify large deletions from gaps in Hi-C coverage.

    Deletions remove chromatin, so the deleted bins have near-zero contact
    counts. The flanking bins on each side of the deletion often show
    elevated contacts as they are brought into proximity.

    Parameters
    ----------
    clr : cooler.Cooler
    region : str
    zero_coverage_threshold : float   fraction of median to call a bin "empty"
    min_gap_bins : int   minimum consecutive empty bins to call a deletion
    flank_enrichment_threshold : float   minimum flanking O/E enrichment

    Returns
    -------
    pd.DataFrame
        chrom, start, end, gap_bins, flank_enrichment, sv_type
    """
    chrom, start, end = parse_coordinates(region)
    matrix = clr.matrix(balance=False).fetch(f"{chrom}:{start}-{end}").astype(float)
    matrix = np.nan_to_num(matrix, nan=0.0)

    coverage = matrix.sum(axis=1)
    median_cov = np.median(coverage[coverage > 0]) if (coverage > 0).any() else 1.0
    empty = coverage < zero_coverage_threshold * median_cov

    resolution = clr.binsize
    results = []
    i = 0
    while i < len(empty):
        if empty[i]:
            j = i
            while j < len(empty) and empty[j]:
                j += 1
            gap_len = j - i
            if gap_len >= min_gap_bins:
                # Check flank enrichment: contacts between left and right flanks
                left_bins = list(range(max(0, i - 5), i))
                right_bins = list(range(j, min(len(empty), j + 5)))
                if left_bins and right_bins:
                    flank_contacts = np.mean([matrix[lb, rb] for lb in left_bins for rb in right_bins])
                    background_contacts = np.mean([matrix[lb, lb + 1] for lb in left_bins[:-1]] or [1.0])
                    enrichment = flank_contacts / max(background_contacts, 1e-10)
                else:
                    enrichment = 0.0

                if enrichment >= flank_enrichment_threshold:
                    results.append({
                        'chrom': chrom,
                        'start': int(start + i * resolution),
                        'end': int(start + j * resolution),
                        'gap_bins': gap_len,
                        'flank_enrichment': float(enrichment),
                        'sv_type': 'deletion',
                    })
            i = j
        else:
            i += 1

    return pd.DataFrame(results) if results else pd.DataFrame(columns=[
        'chrom', 'start', 'end', 'gap_bins', 'flank_enrichment', 'sv_type',
    ])


# ---------------------------------------------------------------------------
# Hi-C based: inversion detection
# ---------------------------------------------------------------------------

def detect_inversions_hic(
    clr: cooler.Cooler,
    region: str,
    window_bins: int = 20,
    di_flip_threshold: float = 2.0,
) -> pd.DataFrame:
    """Detect putative inversions from directionality index sign flips.

    An inversion causes the directionality index (DI) to show an unexpected
    sign pattern: a region that normally biases contacts downstream now
    biases upstream, creating a mirror-image TAD signature.

    Parameters
    ----------
    clr : cooler.Cooler
    region : str
    window_bins : int   DI window size in bins
    di_flip_threshold : float   |DI| threshold to call a flip

    Returns
    -------
    pd.DataFrame
        chrom, start, end, di_left, di_right, flip_score, sv_type
    """
    from .tad_boundaries import compute_directionality_index

    chrom, start, end = parse_coordinates(region)
    resolution = clr.binsize
    window_bp = window_bins * resolution

    di_df = compute_directionality_index(clr, region, window_bp=window_bp)
    di = di_df['DI'].values
    n = len(di)

    results = []
    for i in range(window_bins, n - window_bins):
        left_di = np.nanmean(di[max(0, i - window_bins):i])
        right_di = np.nanmean(di[i:min(n, i + window_bins)])

        # Classic DI boundary: left_di > 0, right_di < 0 (normal TAD boundary)
        # Inversion signature: the entire region has wrong-sign DI
        # Detect: large positive DI on both sides (or both negative) — anomalous
        if (left_di > di_flip_threshold and right_di > di_flip_threshold):
            flip_score = min(left_di, right_di)
            results.append({
                'chrom': chrom,
                'start': int(start + (i - window_bins // 2) * resolution),
                'end': int(start + (i + window_bins // 2) * resolution),
                'di_left': float(left_di),
                'di_right': float(right_di),
                'flip_score': float(flip_score),
                'sv_type': 'inversion_candidate',
            })

    return pd.DataFrame(results) if results else pd.DataFrame(columns=[
        'chrom', 'start', 'end', 'di_left', 'di_right', 'flip_score', 'sv_type',
    ])


# ---------------------------------------------------------------------------
# BAM-based SV detection
# ---------------------------------------------------------------------------

def detect_sv_from_bam(
    bam_path: str,
    region: str,
    min_mapping_quality: int = 20,
    min_support_reads: int = 3,
    max_normal_insert: int = 1000,
    min_split_clip: int = 30,
) -> pd.DataFrame:
    """Detect SVs from discordant and split reads in a BAM file.

    Requires pysam. Two signals are used:
    - Discordant pairs: unexpected orientation (inv/dup) or extreme insert size (del/ins)
    - Split reads: large soft-clips indicating a breakpoint

    Parameters
    ----------
    bam_path : str
    region : str
    min_mapping_quality : int
    min_support_reads : int   minimum read support to report a breakpoint
    max_normal_insert : int   insert sizes above this are discordant
    min_split_clip : int   minimum soft-clip length to count as split

    Returns
    -------
    pd.DataFrame
        chrom, pos, sv_type, support_reads, evidence
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required for BAM-based SV detection. Install with: conda install pysam -c bioconda")

    chrom, start, end = parse_coordinates(region)
    bam = pysam.AlignmentFile(bam_path, "rb")

    discordant_positions: dict = {}
    split_positions: dict = {}

    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.mapping_quality < min_mapping_quality:
            continue
        if read.is_duplicate or read.is_secondary or read.is_supplementary:
            continue

        pos = read.reference_start

        # Discordant pairs
        if read.is_paired and not read.mate_is_unmapped:
            insert = abs(read.template_length)
            if insert > max_normal_insert:
                bin_pos = (pos // 1000) * 1000
                discordant_positions[bin_pos] = discordant_positions.get(bin_pos, 0) + 1
            if read.is_reverse == read.mate_is_reverse:
                bin_pos = (pos // 1000) * 1000
                discordant_positions[bin_pos] = discordant_positions.get(bin_pos, 0) + 1

        # Split reads (large soft clips)
        if read.cigartuples:
            for op, length in read.cigartuples:
                if op == 4 and length >= min_split_clip:  # soft clip
                    bin_pos = (pos // 1000) * 1000
                    split_positions[bin_pos] = split_positions.get(bin_pos, 0) + 1

    bam.close()
    results = []

    for pos, count in discordant_positions.items():
        if count >= min_support_reads:
            results.append({
                'chrom': chrom, 'pos': pos,
                'sv_type': 'discordant_cluster',
                'support_reads': count, 'evidence': 'discordant_pairs',
            })

    for pos, count in split_positions.items():
        if count >= min_support_reads:
            results.append({
                'chrom': chrom, 'pos': pos,
                'sv_type': 'split_read_cluster',
                'support_reads': count, 'evidence': 'split_reads',
            })

    return pd.DataFrame(results) if results else pd.DataFrame(columns=[
        'chrom', 'pos', 'sv_type', 'support_reads', 'evidence',
    ])

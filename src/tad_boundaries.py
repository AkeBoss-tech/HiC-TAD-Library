import cooler
import cooltools
import bioframe
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
from typing import Optional, List, Tuple


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def parse_coordinates(coordinates: str) -> Tuple[str, int, int]:
    """Parse 'chr12:26,000,000-28,000,000' into (chrom, start_bp, end_bp)."""
    chrom, rest = coordinates.split(':')
    start_str, end_str = rest.split('-')
    return chrom, int(start_str.replace(',', '')), int(end_str.replace(',', ''))


def make_view_df(chrom: str, start: int, end: int) -> pd.DataFrame:
    """Single-row view DataFrame for cooltools functions."""
    return pd.DataFrame({
        'chrom': [chrom], 'start': [start], 'end': [end], 'name': [chrom]
    })


# ---------------------------------------------------------------------------
# 1. Directionality Index
# ---------------------------------------------------------------------------

def compute_directionality_index(
    clr: cooler.Cooler,
    coordinates: str,
    window_bp: int = 500_000,
) -> pd.DataFrame:
    """
    Calculate the Directionality Index (Dixon et al.) for each bin in a region.

    DI_i = sign(B-A) * ((A-E)^2 + (B-E)^2) / E
    where A = upstream contacts, B = downstream contacts, E = (A+B)/2.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object at the desired resolution.
    coordinates : str
        Region string, e.g. "chr12:26,000,000-28,000,000".
    window_bp : int
        Size of the upstream/downstream window in base pairs.

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, upstream_sum, downstream_sum, DI.
    """
    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    resolution = clr.binsize
    window_bins = window_bp // resolution

    # Extend fetch region by window size for edge bins
    chrom_len = clr.chromsizes[chrom]
    fetch_start = max(0, start_bp - window_bp)
    fetch_end = min(chrom_len, end_bp + window_bp)
    fetch_coords = f"{chrom}:{fetch_start}-{fetch_end}"

    matrix = clr.matrix(balance=True).fetch(fetch_coords)
    n = matrix.shape[0]

    # Offset of the target region within the fetched matrix
    offset = (start_bp - fetch_start) // resolution
    n_target = (end_bp - start_bp) // resolution

    upstream = np.full(n_target, np.nan)
    downstream = np.full(n_target, np.nan)
    di_values = np.full(n_target, np.nan)

    for idx in range(n_target):
        i = offset + idx
        lo = max(0, i - window_bins)
        hi = min(n, i + window_bins + 1)

        A = np.nansum(matrix[i, lo:i])
        B = np.nansum(matrix[i, i + 1:hi])
        E = (A + B) / 2.0

        upstream[idx] = A
        downstream[idx] = B

        if E == 0:
            di_values[idx] = 0.0
        else:
            di_values[idx] = np.sign(B - A) * ((A - E) ** 2 + (B - E) ** 2) / E

    bins = clr.bins().fetch(f"{chrom}:{start_bp}-{end_bp}")
    result = bins[['chrom', 'start', 'end']].copy()
    result = result.iloc[:n_target].reset_index(drop=True)
    result['upstream_sum'] = upstream
    result['downstream_sum'] = downstream
    result['DI'] = di_values
    return result


# ---------------------------------------------------------------------------
# 2. Boundary prominence scoring
# ---------------------------------------------------------------------------

def _li_threshold(values: np.ndarray, max_iter: int = 1000, tol: float = 1e-6) -> float:
    """Li's iterative minimum cross-entropy thresholding (no scikit-image needed)."""
    vals = values[np.isfinite(values)]
    if len(vals) == 0:
        return 0.0
    t = np.mean(vals)
    for _ in range(max_iter):
        fg = vals[vals > t]
        bg = vals[vals <= t]
        if len(fg) == 0 or len(bg) == 0:
            break
        mfg = np.mean(fg)
        mbg = np.mean(bg)
        if mfg == mbg:
            break
        t_new = (mfg - mbg) / (np.log(mfg) - np.log(mbg)) if mfg > 0 and mbg > 0 else t
        if abs(t_new - t) < tol:
            break
        t = t_new
    return t


def score_boundary_prominence(
    insulation_table: pd.DataFrame,
    window_bp: int,
    prominence_thresholds: Tuple[float, float] = (0.2, 0.5),
    min_distance_bins: int = 3,
    use_auto_threshold: bool = False,
) -> pd.DataFrame:
    """
    Score boundaries by topographic prominence of insulation score minima.

    Parameters
    ----------
    insulation_table : pd.DataFrame
        Output of cooltools.insulation() containing the
        log2_insulation_score_{window_bp} column.
    window_bp : int
        The insulation window size to analyse.
    prominence_thresholds : tuple of (float, float)
        (weak_min, strong_min). Default (0.2, 0.5).
    min_distance_bins : int
        Minimum separation between detected minima.
    use_auto_threshold : bool
        If True, use Li's method to determine the strong/weak threshold
        automatically from the prominence distribution.

    Returns
    -------
    pd.DataFrame
        Original table augmented with columns:
        prominence, boundary_class ("strong", "weak", "sub_threshold", or None).
    """
    ins_col = f'log2_insulation_score_{window_bp}'
    result = insulation_table.copy()
    result['prominence'] = np.nan
    result['boundary_class'] = None

    weak_min, strong_min = prominence_thresholds

    for chrom_name, grp in result.groupby('chrom'):
        ins = grp[ins_col].values.copy()
        # Replace NaN with 0 for peak finding
        mask_nan = np.isnan(ins)
        ins[mask_nan] = 0.0

        # Find minima = peaks of negated signal
        peaks, props = find_peaks(-ins, distance=min_distance_bins, prominence=0)

        if len(peaks) == 0:
            continue

        prominences = props['prominences']

        if use_auto_threshold:
            auto_t = _li_threshold(prominences)
            # Use auto threshold as the weak/strong split
            weak_min = auto_t * 0.5
            strong_min = auto_t

        idx = grp.index[peaks]
        result.loc[idx, 'prominence'] = prominences

        for peak_pos, prom in zip(idx, prominences):
            if prom >= strong_min:
                result.loc[peak_pos, 'boundary_class'] = 'strong'
            elif prom >= weak_min:
                result.loc[peak_pos, 'boundary_class'] = 'weak'
            else:
                result.loc[peak_pos, 'boundary_class'] = 'sub_threshold'

    return result


# ---------------------------------------------------------------------------
# 3. TAD interval calling
# ---------------------------------------------------------------------------

_CLASS_RANK = {'sub_threshold': 0, 'weak': 1, 'strong': 2}


def call_tad_intervals(
    boundary_df: pd.DataFrame,
    clr: cooler.Cooler,
    boundary_class_min: str = 'weak',
    max_tad_length_bp: int = 3_000_000,
    filter_assembly_gaps: bool = False,
    assembly: str = 'mm10',
) -> pd.DataFrame:
    """
    Convert 1D boundary positions to a table of TAD intervals.

    Parameters
    ----------
    boundary_df : pd.DataFrame
        Output of score_boundary_prominence.
    clr : cooler.Cooler
        Cooler object (for bin info and chromsizes).
    boundary_class_min : str
        Minimum boundary class to use as delimiter ("sub_threshold", "weak", "strong").
    max_tad_length_bp : int
        Maximum TAD length. Intervals exceeding this are discarded.
    filter_assembly_gaps : bool
        If True, discard intervals overlapping centromeres/gaps.
    assembly : str
        Genome assembly for gap annotation.

    Returns
    -------
    pd.DataFrame
        BED-like with columns: chrom, start, end, tad_id, length_bp, n_bins,
        left_boundary_class, right_boundary_class.
    """
    min_rank = _CLASS_RANK.get(boundary_class_min, 1)

    # Filter to qualifying boundaries
    boundaries = boundary_df[
        boundary_df['boundary_class'].map(lambda c: _CLASS_RANK.get(c, -1) >= min_rank)
    ].copy()

    resolution = clr.binsize
    tads = []
    tad_id = 0

    for chrom_name, chrom_bounds in boundaries.groupby('chrom'):
        chrom_len = clr.chromsizes[chrom_name]
        sorted_bounds = chrom_bounds.sort_values('start')

        positions = sorted_bounds[['start', 'end', 'boundary_class']].values.tolist()

        # First TAD: chromosome start → first boundary
        if len(positions) > 0:
            first_bnd_start, first_bnd_end, first_cls = positions[0]
            if first_bnd_start > 0:
                tads.append({
                    'chrom': chrom_name,
                    'start': 0,
                    'end': int(first_bnd_start),
                    'tad_id': tad_id,
                    'left_boundary_class': None,
                    'right_boundary_class': first_cls,
                })
                tad_id += 1

        # Interior TADs: between consecutive boundaries
        for j in range(len(positions) - 1):
            left_start, left_end, left_cls = positions[j]
            right_start, right_end, right_cls = positions[j + 1]
            tads.append({
                'chrom': chrom_name,
                'start': int(left_end),
                'end': int(right_start),
                'tad_id': tad_id,
                'left_boundary_class': left_cls,
                'right_boundary_class': right_cls,
            })
            tad_id += 1

        # Last TAD: last boundary → chromosome end
        if len(positions) > 0:
            last_start, last_end, last_cls = positions[-1]
            if last_end < chrom_len:
                tads.append({
                    'chrom': chrom_name,
                    'start': int(last_end),
                    'end': int(chrom_len),
                    'tad_id': tad_id,
                    'left_boundary_class': last_cls,
                    'right_boundary_class': None,
                })
                tad_id += 1

    if not tads:
        return pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'tad_id', 'length_bp', 'n_bins',
            'left_boundary_class', 'right_boundary_class'
        ])

    tad_df = pd.DataFrame(tads)
    tad_df['length_bp'] = tad_df['end'] - tad_df['start']
    tad_df['n_bins'] = tad_df['length_bp'] // resolution

    # Filter by maximum length
    tad_df = tad_df[tad_df['length_bp'] <= max_tad_length_bp].copy()

    # Filter by assembly gaps
    if filter_assembly_gaps:
        try:
            gaps = bioframe.fetch_centromeres(assembly)
            overlaps = bioframe.overlap(tad_df, gaps, how='inner')
            tad_df = tad_df[~tad_df.index.isin(overlaps.index)]
        except Exception:
            print(f"  Warning: could not fetch gap annotations for {assembly}, skipping gap filter.")

    tad_df = tad_df.reset_index(drop=True)
    tad_df['tad_id'] = range(len(tad_df))
    return tad_df


# ---------------------------------------------------------------------------
# 4. Multi-scale insulation
# ---------------------------------------------------------------------------

def multiscale_insulation(
    clr: cooler.Cooler,
    coordinates: str,
    window_sizes_bp: Optional[List[int]] = None,
    n_windows: int = 8,
    ignore_diags: int = 2,
) -> Tuple[pd.DataFrame, List[int]]:
    """
    Run insulation analysis at multiple window sizes (log-spaced by default).

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object at the desired resolution.
    coordinates : str
        Region string.
    window_sizes_bp : list of int or None
        Explicit window sizes. If None, auto-generates log-spaced sizes.
    n_windows : int
        Number of log-spaced windows when auto-generating.
    ignore_diags : int
        Diagonals to ignore in the insulation calculation.

    Returns
    -------
    (pd.DataFrame, list of int)
        The insulation table filtered to the target region, and the list of
        window sizes used.
    """
    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    resolution = clr.binsize

    if window_sizes_bp is None:
        raw = np.geomspace(25_000, 500_000, n_windows)
        window_sizes_bp = sorted(set(
            int(round(w / resolution) * resolution) for w in raw
        ))
        # Ensure all windows are at least 2 * resolution
        window_sizes_bp = [w for w in window_sizes_bp if w >= 2 * resolution]

    # Context buffer = largest window
    max_window = max(window_sizes_bp)
    chrom_len = clr.chromsizes[chrom]
    fetch_start = max(0, start_bp - max_window)
    fetch_end = min(chrom_len, end_bp + max_window)
    view_df = make_view_df(chrom, fetch_start, fetch_end)

    ins_table = cooltools.insulation(
        clr, window_sizes_bp, ignore_diags=ignore_diags, view_df=view_df
    )

    # Filter to target region
    region_ins = ins_table[
        (ins_table['chrom'] == chrom) &
        (ins_table['start'] >= start_bp) &
        (ins_table['end'] <= end_bp)
    ].copy().reset_index(drop=True)

    return region_ins, window_sizes_bp


def classify_boundary_persistence(
    multiscale_table: pd.DataFrame,
    window_sizes_bp: List[int],
    prominence_threshold: float = 0.2,
    min_distance_bins: int = 3,
) -> pd.DataFrame:
    """
    Classify boundaries by persistence across insulation scales.

    Parameters
    ----------
    multiscale_table : pd.DataFrame
        Output of multiscale_insulation.
    window_sizes_bp : list of int
        Window sizes used.
    prominence_threshold : float
        Minimum prominence to count as a boundary at each scale.
    min_distance_bins : int
        Minimum peak separation.

    Returns
    -------
    pd.DataFrame
        Augmented with: n_scales_boundary, fraction_scales, boundary_type
        ("constitutive", "nested_sub", "nested_meta", or None),
        min_scale_bp, max_scale_bp.
    """
    result = multiscale_table[['chrom', 'start', 'end']].copy()
    n_scales = len(window_sizes_bp)
    median_window = np.median(window_sizes_bp)

    # Track which scales detect each bin as a boundary
    boundary_at_scale = np.zeros((len(result), n_scales), dtype=bool)
    scale_detected = [[] for _ in range(len(result))]

    for si, w in enumerate(window_sizes_bp):
        ins_col = f'log2_insulation_score_{w}'
        if ins_col not in multiscale_table.columns:
            continue

        scored = score_boundary_prominence(
            multiscale_table, w,
            prominence_thresholds=(prominence_threshold, prominence_threshold * 2),
            min_distance_bins=min_distance_bins,
        )
        is_boundary = scored['boundary_class'].notna()
        boundary_at_scale[:len(is_boundary), si] = is_boundary.values[:len(result)]

        for idx in np.where(is_boundary.values[:len(result)])[0]:
            scale_detected[idx].append(w)

    result['n_scales_boundary'] = boundary_at_scale.sum(axis=1)
    result['fraction_scales'] = result['n_scales_boundary'] / n_scales

    types = []
    min_scales = []
    max_scales = []
    for idx in range(len(result)):
        detected = scale_detected[idx]
        if not detected:
            types.append(None)
            min_scales.append(np.nan)
            max_scales.append(np.nan)
        else:
            min_scales.append(min(detected))
            max_scales.append(max(detected))
            frac = result.iloc[idx]['fraction_scales']
            if frac >= 0.75:
                types.append('constitutive')
            elif all(s < median_window for s in detected):
                types.append('nested_sub')
            elif all(s >= median_window for s in detected):
                types.append('nested_meta')
            else:
                types.append('mixed')

    result['boundary_type'] = types
    result['min_scale_bp'] = min_scales
    result['max_scale_bp'] = max_scales
    return result


# ---------------------------------------------------------------------------
# 5. Boundary pileup
# ---------------------------------------------------------------------------

def boundary_pileup(
    clr: cooler.Cooler,
    boundary_df: pd.DataFrame,
    flank_bins: int = 40,
    boundary_class_min: str = 'weak',
    normalize_by_expected: bool = True,
) -> Tuple[np.ndarray, int]:
    """
    Aggregate 2D contact sub-matrices centered on boundaries.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object.
    boundary_df : pd.DataFrame
        Output of score_boundary_prominence.
    flank_bins : int
        Number of bins on each side of the boundary center.
    boundary_class_min : str
        Minimum boundary class to include.
    normalize_by_expected : bool
        If True, divide each pixel by distance-dependent expected.

    Returns
    -------
    (np.ndarray, int)
        The averaged pileup matrix of shape (2*flank_bins+1, 2*flank_bins+1),
        and the number of boundaries that were stacked.
    """
    min_rank = _CLASS_RANK.get(boundary_class_min, 1)
    boundaries = boundary_df[
        boundary_df['boundary_class'].map(lambda c: _CLASS_RANK.get(c, -1) >= min_rank)
    ]

    resolution = clr.binsize
    snippet_size = 2 * flank_bins + 1

    # Precompute expected per chromosome if normalizing
    expected_dict = {}
    if normalize_by_expected:
        for chrom_name in boundaries['chrom'].unique():
            chrom_len = clr.chromsizes[chrom_name]
            view_df = make_view_df(chrom_name, 0, chrom_len)
            try:
                exp = cooltools.expected_cis(clr, view_df=view_df)
                exp_values = exp['balanced.avg'].values
                expected_dict[chrom_name] = exp_values
            except Exception:
                expected_dict[chrom_name] = None

    stack = []
    for _, row in boundaries.iterrows():
        chrom_name = row['chrom']
        mid_bp = (row['start'] + row['end']) // 2
        center_bin = mid_bp // resolution

        # Bin range within the chromosome
        chrom_n_bins = clr.chromsizes[chrom_name] // resolution
        lo = center_bin - flank_bins
        hi = center_bin + flank_bins + 1
        if lo < 0 or hi > chrom_n_bins:
            continue

        start_bp = lo * resolution
        end_bp_fetch = hi * resolution
        snippet = clr.matrix(balance=True).fetch(
            f"{chrom_name}:{start_bp}-{end_bp_fetch}"
        )

        if snippet.shape != (snippet_size, snippet_size):
            continue

        if normalize_by_expected and expected_dict.get(chrom_name) is not None:
            exp_vals = expected_dict[chrom_name]
            exp_matrix = np.zeros_like(snippet)
            for i in range(snippet_size):
                for j in range(snippet_size):
                    d = abs(j - i)
                    if d < len(exp_vals):
                        exp_matrix[i, j] = exp_vals[d]
            with np.errstate(divide='ignore', invalid='ignore'):
                snippet = np.where(exp_matrix > 0, snippet / exp_matrix, np.nan)

        stack.append(snippet)

    if not stack:
        return np.full((snippet_size, snippet_size), np.nan), 0

    pileup = np.nanmean(np.array(stack), axis=0)
    return pileup, len(stack)

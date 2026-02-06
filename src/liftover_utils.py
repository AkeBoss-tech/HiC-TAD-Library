"""
Coordinate liftover utilities for cross-species TAD boundary comparison.

Converts mm10 (mouse) TAD boundary coordinates to hg38 (human) using pyliftover.
Requires: pip install pyliftover
"""

import numpy as np
import pandas as pd
from typing import Optional, Tuple


def liftover_boundaries(
    boundary_df: pd.DataFrame,
    source_assembly: str = 'mm10',
    target_assembly: str = 'hg38',
) -> pd.DataFrame:
    """
    Lift over boundary coordinates from one genome assembly to another.

    Parameters
    ----------
    boundary_df : pd.DataFrame
        Must have columns: chrom, start, end.
        Optionally: boundary_class, prominence.
    source_assembly : str
        Source genome assembly (default: mm10).
    target_assembly : str
        Target genome assembly (default: hg38).

    Returns
    -------
    pd.DataFrame
        One row per input boundary with source coords, target coords,
        and liftover_success flag.
    """
    from pyliftover import LiftOver

    lo = LiftOver(source_assembly, target_assembly)

    results = []
    for _, row in boundary_df.iterrows():
        chrom = row['chrom']
        mid = (int(row['start']) + int(row['end'])) // 2

        converted = lo.convert_coordinate(chrom, mid)

        base = {
            'source_chrom': chrom,
            'source_start': int(row['start']),
            'source_end': int(row['end']),
            'source_mid': mid,
            'boundary_class': row.get('boundary_class'),
            'prominence': row.get('prominence'),
        }

        if converted and len(converted) > 0:
            target_chrom, target_pos, target_strand, _ = converted[0]
            base.update({
                'target_chrom': target_chrom,
                'target_pos': int(target_pos),
                'target_strand': target_strand,
                'liftover_success': True,
            })
        else:
            base.update({
                'target_chrom': None,
                'target_pos': None,
                'target_strand': None,
                'liftover_success': False,
            })

        results.append(base)

    return pd.DataFrame(results)


def liftover_region(
    coordinates: str,
    source_assembly: str = 'mm10',
    target_assembly: str = 'hg38',
) -> Optional[str]:
    """
    Lift over a region string (e.g. 'chr12:26,000,000-28,000,000') to
    another assembly.  Returns a target region of the same size centered
    on the lifted-over midpoint, or None if the midpoint cannot be mapped.
    """
    from pyliftover import LiftOver

    lo = LiftOver(source_assembly, target_assembly)

    chrom, rest = coordinates.split(':')
    start_str, end_str = rest.split('-')
    start = int(start_str.replace(',', ''))
    end = int(end_str.replace(',', ''))
    region_size = end - start
    mid = (start + end) // 2

    converted = lo.convert_coordinate(chrom, mid)

    if converted and len(converted) > 0:
        target_chrom, target_pos, _, _ = converted[0]
        target_start = max(0, int(target_pos) - region_size // 2)
        target_end = int(target_pos) + region_size // 2
        return f"{target_chrom}:{target_start}-{target_end}"

    return None

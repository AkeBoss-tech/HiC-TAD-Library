"""SNV and small indel calling from BAM files.

Implements a lightweight pileup-based variant caller using pysam.
This is not a replacement for GATK/DeepVariant but provides a simple
framework for exploratory analysis or validating known variant positions.

Algorithm:
    1. Pileup reads over a region.
    2. At each position, count A/T/G/C/del observations filtered by
       base quality and mapping quality.
    3. Call a variant when:
       - At least min_depth reads cover the position.
       - The most common alt allele frequency (VAF) exceeds min_vaf.
       - The variant is not the reference allele.
    4. Annotate each call with depth, VAF, strand bias, and a simple
       quality score.
"""

import numpy as np
import pandas as pd
from typing import Optional

from .tad_boundaries import parse_coordinates

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False


_BASES = ('A', 'C', 'G', 'T')


# ---------------------------------------------------------------------------
# Core pileup caller
# ---------------------------------------------------------------------------

def call_snvs(
    bam_path: str,
    region: str,
    ref_fasta: Optional[str] = None,
    min_depth: int = 10,
    min_vaf: float = 0.20,
    min_base_quality: int = 20,
    min_mapping_quality: int = 30,
    strand_bias_threshold: float = 0.1,
) -> pd.DataFrame:
    """Call SNVs and small indels in a region from a BAM file.

    Parameters
    ----------
    bam_path : str
    region : str   e.g. 'chr12:26,000,000-28,000,000'
    ref_fasta : str, optional   path to indexed FASTA (enables ref-allele annotation)
    min_depth : int   minimum read depth
    min_vaf : float   minimum variant allele frequency (0–1)
    min_base_quality : int
    min_mapping_quality : int
    strand_bias_threshold : float   maximum |forward_frac - 0.5| to flag strand bias

    Returns
    -------
    pd.DataFrame
        chrom, pos, ref, alt, depth, alt_count, vaf, strand_bias, qual, variant_type
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required for variant calling. Install with: conda install pysam -c bioconda")

    chrom, start, end = parse_coordinates(region)

    bam = pysam.AlignmentFile(bam_path, "rb")
    fasta = pysam.FastaFile(ref_fasta) if ref_fasta else None

    calls = []

    for pileup_col in bam.pileup(
        chrom, start, end,
        min_base_quality=min_base_quality,
        min_mapping_quality=min_mapping_quality,
        truncate=True,
        ignore_overlaps=True,
    ):
        pos = pileup_col.reference_pos
        if pos < start or pos >= end:
            continue

        base_counts = {b: 0 for b in _BASES}
        del_count = 0
        ins_count = 0
        fwd_counts = {b: 0 for b in _BASES}

        for pread in pileup_col.pileups:
            read = pread.alignment
            if pread.is_del:
                del_count += 1
                continue
            if pread.is_refskip:
                continue
            qpos = pread.query_position
            if qpos is None:
                continue
            base = read.query_sequence[qpos].upper()
            if base in _BASES:
                base_counts[base] += 1
                if not read.is_reverse:
                    fwd_counts[base] += 1

        depth = sum(base_counts.values()) + del_count
        if depth < min_depth:
            continue

        # Reference allele
        ref_base = 'N'
        if fasta:
            try:
                ref_base = fasta.fetch(chrom, pos, pos + 1).upper()
            except Exception:
                pass

        # Find most common alt allele
        alt_counts_sorted = sorted(
            [(b, c) for b, c in base_counts.items() if b != ref_base],
            key=lambda x: -x[1],
        )
        if not alt_counts_sorted:
            continue

        alt_base, alt_count = alt_counts_sorted[0]
        vaf = alt_count / depth
        if vaf < min_vaf:
            continue

        # Strand bias
        total_fwd = fwd_counts[alt_base]
        strand_bias = abs(total_fwd / max(alt_count, 1) - 0.5)
        sb_flag = strand_bias > (0.5 - strand_bias_threshold)

        # Simple quality score: depth * vaf * (1 - strand_bias)
        qual = int(depth * vaf * (1.0 - strand_bias))

        # Variant type
        if ref_base == 'N':
            vtype = 'snv'
        elif len(ref_base) == len(alt_base):
            vtype = 'snv'
        else:
            vtype = 'indel'

        calls.append({
            'chrom': chrom,
            'pos': int(pos),
            'ref': ref_base,
            'alt': alt_base,
            'depth': int(depth),
            'alt_count': int(alt_count),
            'vaf': float(vaf),
            'strand_bias': float(strand_bias),
            'strand_bias_flag': sb_flag,
            'qual': qual,
            'variant_type': vtype,
        })

    if fasta:
        fasta.close()
    bam.close()

    if not calls:
        return pd.DataFrame(columns=[
            'chrom', 'pos', 'ref', 'alt', 'depth', 'alt_count',
            'vaf', 'strand_bias', 'strand_bias_flag', 'qual', 'variant_type',
        ])

    df = pd.DataFrame(calls).sort_values('pos').reset_index(drop=True)
    return df


# ---------------------------------------------------------------------------
# VCF parser (for downstream phasing / annotation)
# ---------------------------------------------------------------------------

def parse_vcf(
    vcf_path: str,
    region: Optional[str] = None,
    min_qual: float = 0.0,
    pass_only: bool = True,
) -> pd.DataFrame:
    """Parse a VCF file into a DataFrame.

    Parameters
    ----------
    vcf_path : str
    region : str, optional   restrict to region
    min_qual : float   minimum QUAL value
    pass_only : bool   only include PASS variants

    Returns
    -------
    pd.DataFrame
        chrom, pos, ref, alt, qual, filter, genotype, is_het
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required for VCF parsing.")

    vcf = pysam.VariantFile(vcf_path)

    chrom_filter, start_filter, end_filter = None, None, None
    if region:
        chrom_filter, start_filter, end_filter = parse_coordinates(region)

    records = []
    fetch_kwargs = {}
    if chrom_filter:
        fetch_kwargs = {'contig': chrom_filter, 'start': start_filter, 'stop': end_filter}

    for rec in vcf.fetch(**fetch_kwargs):
        if pass_only and rec.filter.keys() and 'PASS' not in rec.filter:
            continue
        if rec.qual is not None and rec.qual < min_qual:
            continue

        for sample_name, sample in rec.samples.items():
            gt = sample.get('GT', (None, None))
            is_het = len(set(g for g in gt if g is not None)) > 1 if gt else False
            records.append({
                'chrom': rec.chrom,
                'pos': rec.pos,
                'ref': rec.ref,
                'alt': ','.join(str(a) for a in rec.alts) if rec.alts else '.',
                'qual': float(rec.qual) if rec.qual is not None else np.nan,
                'filter': ','.join(rec.filter.keys()),
                'sample': sample_name,
                'genotype': '/'.join(str(g) for g in gt) if gt else './.',
                'is_het': is_het,
            })
            break  # first sample only

    vcf.close()
    return pd.DataFrame(records) if records else pd.DataFrame(columns=[
        'chrom', 'pos', 'ref', 'alt', 'qual', 'filter', 'sample', 'genotype', 'is_het',
    ])


# ---------------------------------------------------------------------------
# Transition/transversion annotation
# ---------------------------------------------------------------------------

_TRANSITIONS = {('A', 'G'), ('G', 'A'), ('C', 'T'), ('T', 'C')}


def annotate_ti_tv(df: pd.DataFrame) -> pd.DataFrame:
    """Add Ti/Tv classification column to a SNV DataFrame."""
    result = df.copy()

    def _classify(row):
        r, a = str(row.get('ref', '')).upper(), str(row.get('alt', '')).upper()
        if len(r) != 1 or len(a) != 1:
            return 'indel'
        if (r, a) in _TRANSITIONS:
            return 'transition'
        if r in _BASES and a in _BASES:
            return 'transversion'
        return 'unknown'

    result['ti_tv'] = result.apply(_classify, axis=1)
    snv_mask = result['ti_tv'].isin(['transition', 'transversion'])
    if snv_mask.any():
        ti = (result.loc[snv_mask, 'ti_tv'] == 'transition').sum()
        tv = (result.loc[snv_mask, 'ti_tv'] == 'transversion').sum()
        result.attrs['ti_tv_ratio'] = ti / tv if tv > 0 else np.nan
    return result

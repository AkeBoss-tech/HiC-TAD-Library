"""Haplotype phasing using Hi-C read pairs.

Hi-C read pairs physically connect distant genomic loci within the same
nucleus. When a read pair spans two heterozygous SNP positions, the
alleles observed at each SNP must come from the same chromosome (cis
configuration), providing a phasing constraint.

Algorithm (simplified Hi-C phasing, as in HapCUT2-HiC):
    1. Load heterozygous SNPs from a VCF file.
    2. Scan a BAM file for read pairs in which at least one read overlaps
       a het SNP and the paired read overlaps another.
    3. Build a phasing graph: nodes = SNP alleles, edges = co-occurrence
       in the same read pair with a weight = number of supporting pairs.
    4. Resolve the graph into haplotype blocks by greedy bipartition of
       connected components.
    5. Report phased haplotype blocks and individual SNP assignments.

Note: This module requires pysam.  For production phasing use
      HapCUT2, WhatsHap, or 3D-DNA; this implementation is a
      self-contained educational version.
"""

import numpy as np
import pandas as pd
from collections import defaultdict
from typing import Dict, List, Optional, Tuple

from .tad_boundaries import parse_coordinates

try:
    import pysam
    _PYSAM_AVAILABLE = True
except ImportError:
    _PYSAM_AVAILABLE = False


# ---------------------------------------------------------------------------
# SNP loading
# ---------------------------------------------------------------------------

def load_het_snps(
    vcf_path: str,
    region: Optional[str] = None,
) -> pd.DataFrame:
    """Load heterozygous SNPs from a VCF file.

    Parameters
    ----------
    vcf_path : str
    region : str, optional   restrict to region

    Returns
    -------
    pd.DataFrame  columns: chrom, pos, ref, alt, snp_id
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required for VCF parsing.")

    vcf = pysam.VariantFile(vcf_path)
    fetch_kwargs = {}
    if region:
        chrom, start, end = parse_coordinates(region)
        fetch_kwargs = {'contig': chrom, 'start': start, 'stop': end}

    records = []
    for rec in vcf.fetch(**fetch_kwargs):
        if rec.alts is None:
            continue
        for sample_name, sample in rec.samples.items():
            gt = sample.get('GT', (None, None))
            if gt is None:
                continue
            alleles = set(g for g in gt if g is not None)
            if len(alleles) == 2:  # heterozygous
                records.append({
                    'chrom': rec.chrom,
                    'pos': rec.pos,
                    'ref': rec.ref,
                    'alt': rec.alts[0],
                    'snp_id': f"{rec.chrom}:{rec.pos}",
                })
            break  # first sample

    vcf.close()
    df = pd.DataFrame(records)
    if df.empty:
        return pd.DataFrame(columns=['chrom', 'pos', 'ref', 'alt', 'snp_id'])
    return df.sort_values(['chrom', 'pos']).reset_index(drop=True)


# ---------------------------------------------------------------------------
# Read-pair to SNP overlap
# ---------------------------------------------------------------------------

def _base_at_pos(read: 'pysam.AlignedSegment', ref_pos: int) -> Optional[str]:
    """Return the query base at a reference position, or None."""
    for qpos, rpos in read.get_aligned_pairs():
        if rpos == ref_pos and qpos is not None:
            return read.query_sequence[qpos].upper()
    return None


def extract_phasing_constraints(
    bam_path: str,
    snps_df: pd.DataFrame,
    region: str,
    min_mapq: int = 20,
    min_base_quality: int = 20,
) -> pd.DataFrame:
    """Find read pairs that span two het SNPs and extract phasing links.

    Parameters
    ----------
    bam_path : str
    snps_df : pd.DataFrame   output of load_het_snps
    region : str
    min_mapq : int
    min_base_quality : int

    Returns
    -------
    pd.DataFrame
        snp1_id, snp1_allele, snp2_id, snp2_allele, read_pair_id, same_haplotype
    """
    if not _PYSAM_AVAILABLE:
        raise ImportError("pysam is required.")

    chrom, start, end = parse_coordinates(region)

    # Index SNPs by position for fast lookup
    region_snps = snps_df[snps_df['chrom'] == chrom].copy()
    snp_by_pos: Dict[int, dict] = {
        int(row['pos']): row.to_dict() for _, row in region_snps.iterrows()
    }
    snp_positions = sorted(snp_by_pos.keys())

    if not snp_positions:
        return pd.DataFrame(columns=[
            'snp1_id', 'snp1_allele', 'snp2_id', 'snp2_allele',
            'read_pair_id', 'same_haplotype',
        ])

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Group reads by query name (pair identification)
    read_groups: Dict[str, List] = defaultdict(list)
    for read in bam.fetch(chrom, start, end):
        if read.is_unmapped or read.mapping_quality < min_mapq:
            continue
        if read.is_duplicate or read.is_secondary:
            continue
        read_groups[read.query_name].append(read)

    bam.close()

    constraints = []
    for read_name, reads in read_groups.items():
        # Find SNPs covered by each read in the pair
        read_snp_alleles: List[Tuple[str, str]] = []  # (snp_id, allele)

        for read in reads:
            if read.reference_start is None:
                continue
            r_start = read.reference_start
            r_end = read.reference_end or r_start

            for snp_pos in snp_positions:
                if not (r_start <= snp_pos < r_end):
                    continue
                snp_info = snp_by_pos[snp_pos]
                base = _base_at_pos(read, snp_pos)
                if base is None:
                    continue
                if base == snp_info['ref']:
                    allele = 'REF'
                elif base == snp_info['alt']:
                    allele = 'ALT'
                else:
                    continue  # unknown base
                read_snp_alleles.append((snp_info['snp_id'], allele))

        # Pairs of SNPs observed in the same read pair = phasing constraint
        for k1 in range(len(read_snp_alleles)):
            for k2 in range(k1 + 1, len(read_snp_alleles)):
                s1_id, s1_allele = read_snp_alleles[k1]
                s2_id, s2_allele = read_snp_alleles[k2]
                if s1_id == s2_id:
                    continue
                # Same read pair → same haplotype if same allele type (REF-REF or ALT-ALT)
                same_hap = s1_allele == s2_allele
                constraints.append({
                    'snp1_id': s1_id,
                    'snp1_allele': s1_allele,
                    'snp2_id': s2_id,
                    'snp2_allele': s2_allele,
                    'read_pair_id': read_name,
                    'same_haplotype': same_hap,
                })

    return pd.DataFrame(constraints) if constraints else pd.DataFrame(columns=[
        'snp1_id', 'snp1_allele', 'snp2_id', 'snp2_allele',
        'read_pair_id', 'same_haplotype',
    ])


# ---------------------------------------------------------------------------
# Graph-based phasing
# ---------------------------------------------------------------------------

def phase_snps(
    constraints_df: pd.DataFrame,
    snps_df: pd.DataFrame,
    min_support: int = 2,
) -> pd.DataFrame:
    """Phase SNPs into haplotype blocks from pairwise constraints.

    Parameters
    ----------
    constraints_df : pd.DataFrame   output of extract_phasing_constraints
    snps_df : pd.DataFrame          het SNPs
    min_support : int               minimum pair support to trust a link

    Returns
    -------
    pd.DataFrame
        snp_id, chrom, pos, ref, alt, haplotype_block, haplotype (0 or 1)
    """
    if constraints_df.empty:
        result = snps_df[['snp_id', 'chrom', 'pos', 'ref', 'alt']].copy()
        result['haplotype_block'] = -1
        result['haplotype'] = -1
        return result

    # Count concordant (same hap) and discordant links
    link_counts: Dict[Tuple[str, str], Dict[str, int]] = defaultdict(lambda: {'same': 0, 'diff': 0})
    for _, row in constraints_df.iterrows():
        key = tuple(sorted([row['snp1_id'], row['snp2_id']]))
        if row['same_haplotype']:
            link_counts[key]['same'] += 1
        else:
            link_counts[key]['diff'] += 1

    # Build weighted phasing graph
    # Nodes: SNP IDs. Edges: (same/diff, weight)
    snp_ids = list(snps_df['snp_id'].unique())
    snp_idx = {s: i for i, s in enumerate(snp_ids)}
    n = len(snp_ids)

    adj: Dict[int, List[Tuple[int, bool]]] = defaultdict(list)
    for (s1, s2), counts in link_counts.items():
        if s1 not in snp_idx or s2 not in snp_idx:
            continue
        total = counts['same'] + counts['diff']
        if total < min_support:
            continue
        same = counts['same'] > counts['diff']
        i1, i2 = snp_idx[s1], snp_idx[s2]
        adj[i1].append((i2, same))
        adj[i2].append((i1, same))

    # BFS to assign haplotypes within connected components
    haplotype = np.full(n, -1, dtype=int)
    block_id = np.full(n, -1, dtype=int)
    current_block = 0

    for start_node in range(n):
        if haplotype[start_node] != -1:
            continue
        if not adj[start_node]:
            haplotype[start_node] = 0
            block_id[start_node] = current_block
            current_block += 1
            continue

        # BFS
        queue = [start_node]
        haplotype[start_node] = 0
        block_id[start_node] = current_block
        head = 0

        while head < len(queue):
            node = queue[head]
            head += 1
            for neighbour, same_hap in adj[node]:
                assigned_hap = haplotype[node] if same_hap else (1 - haplotype[node])
                if haplotype[neighbour] == -1:
                    haplotype[neighbour] = assigned_hap
                    block_id[neighbour] = current_block
                    queue.append(neighbour)
        current_block += 1

    # Build output
    result = snps_df[['snp_id', 'chrom', 'pos', 'ref', 'alt']].copy().reset_index(drop=True)
    hap_map = {snp_ids[i]: (int(haplotype[i]), int(block_id[i])) for i in range(n)}
    result['haplotype_block'] = result['snp_id'].map(lambda s: hap_map.get(s, (-1, -1))[1])
    result['haplotype'] = result['snp_id'].map(lambda s: hap_map.get(s, (-1, -1))[0])
    return result


# ---------------------------------------------------------------------------
# Convenience wrapper
# ---------------------------------------------------------------------------

def run_hic_phasing(
    bam_path: str,
    vcf_path: str,
    region: str,
    min_mapq: int = 20,
    min_support: int = 2,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run the full Hi-C phasing pipeline for a region.

    Returns
    -------
    (phased_snps_df, constraints_df)
    """
    snps = load_het_snps(vcf_path, region=region)
    if snps.empty:
        print(f"  Phasing: no het SNPs found in {region}")
        return snps, pd.DataFrame()

    constraints = extract_phasing_constraints(bam_path, snps, region, min_mapq=min_mapq)
    phased = phase_snps(constraints, snps, min_support=min_support)
    return phased, constraints

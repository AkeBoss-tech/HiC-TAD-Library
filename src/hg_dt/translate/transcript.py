import pandas as pd
from typing import Dict, List
from Bio.Seq import Seq

def extract_exons(genes_df: pd.DataFrame, transcript_id: str) -> pd.DataFrame:
    """Extract and sort exons for a specific transcript."""
    if genes_df.empty or 'Feature' not in genes_df.columns:
        return pd.DataFrame()

    exons = genes_df[(genes_df['Feature'] == 'exon') & (genes_df['transcript_id'] == transcript_id)]

    # Sort by start coordinate.
    exons = exons.sort_values(by='Start')
    return exons

def get_transcript_ids(genes_df: pd.DataFrame) -> List[str]:
    """Get all unique transcript IDs."""
    if genes_df.empty or 'transcript_id' not in genes_df.columns:
        return []
    return list(genes_df['transcript_id'].dropna().unique())

def get_mut_coord(genomic_coord: int, win_start: int, edit_start: int, edit_end: int, mut_shift: int) -> int:
    """Map a genomic coordinate to the relative coordinate in the mutant sequence window."""
    if genomic_coord < edit_start:
        return genomic_coord - win_start
    elif genomic_coord < edit_end:
        # Inside the edit (e.g. deleted region)
        return edit_start - win_start
    else:
        return genomic_coord - win_start + mut_shift

def predict_isoforms(genes_df: pd.DataFrame, ref_seq: str, mut_seq: str, win_start: int,
                     edit_start: int, edit_end: int, mut_shift: int) -> Dict[str, Dict[str, str]]:
    """
    Predict the mature mRNA(s) following DNA edits.

    Args:
        genes_df: The GTF annotations DataFrame.
        ref_seq: The reference sequence window.
        mut_seq: The mutant sequence window.
        win_start: The genomic start position of the window.
        edit_start: Genomic start position of the edit.
        edit_end: Genomic end position of the edit.
        mut_shift: The net change in sequence length (e.g. -1 for 1bp deletion).

    Returns:
        A dictionary mapping transcript_id to a dict containing 'ref_mrna' and 'mut_mrna'.
    """
    isoforms = {}
    transcripts = get_transcript_ids(genes_df)

    for tid in transcripts:
        exons = extract_exons(genes_df, tid)
        if exons.empty:
            continue

        # Get the strand. We assume all exons of a transcript are on the same strand.
        strand = exons.iloc[0]['Strand']

        ref_mrna_dna = ""
        mut_mrna_dna = ""

        for _, exon in exons.iterrows():
            exon_start = exon['Start']
            exon_end = exon['End']

            # Reference assembly
            rel_start = exon_start - win_start
            rel_end = exon_end - win_start

            if rel_start < 0: rel_start = 0
            if rel_end > len(ref_seq): rel_end = len(ref_seq)

            if rel_start < rel_end:
                ref_mrna_dna += ref_seq[rel_start:rel_end]

            # Mutant assembly
            mut_start = get_mut_coord(exon_start, win_start, edit_start, edit_end, mut_shift)
            mut_end = get_mut_coord(exon_end, win_start, edit_start, edit_end, mut_shift)

            if mut_start < 0: mut_start = 0
            if mut_end > len(mut_seq): mut_end = len(mut_seq)

            if mut_start < mut_end:
                mut_mrna_dna += mut_seq[mut_start:mut_end]

        # Reverse complement if on negative strand
        if strand == '-':
            ref_mrna_dna = str(Seq(ref_mrna_dna).reverse_complement())
            mut_mrna_dna = str(Seq(mut_mrna_dna).reverse_complement())

        # Convert to RNA
        ref_mrna = ref_mrna_dna.upper().replace('T', 'U')
        mut_mrna = mut_mrna_dna.upper().replace('T', 'U')

        isoforms[tid] = {
            'ref_mrna': ref_mrna,
            'mut_mrna': mut_mrna
        }

    return isoforms

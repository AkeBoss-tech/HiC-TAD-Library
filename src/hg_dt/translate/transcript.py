import numpy as np
import pandas as pd
from typing import Dict, List, Optional
from Bio.Seq import Seq


def extract_exons(genes_df: pd.DataFrame, transcript_id: str) -> pd.DataFrame:
    """Extract and sort exons for a specific transcript."""
    if genes_df.empty or 'Feature' not in genes_df.columns:
        return pd.DataFrame()
    exons = genes_df[(genes_df['Feature'] == 'exon') & (genes_df['transcript_id'] == transcript_id)]
    return exons.sort_values(by='Start')


def get_transcript_ids(genes_df: pd.DataFrame) -> List[str]:
    """Get all unique transcript IDs."""
    if genes_df.empty or 'transcript_id' not in genes_df.columns:
        return []
    return list(genes_df['transcript_id'].dropna().unique())


def get_mut_coord(genomic_coord: int, win_start: int, edit_start: int,
                  edit_end: int, mut_shift: int) -> int:
    """Map a genomic coordinate to its relative coordinate in the mutant window."""
    if genomic_coord < edit_start:
        return genomic_coord - win_start
    elif genomic_coord < edit_end:
        return edit_start - win_start
    else:
        return genomic_coord - win_start + mut_shift


def _score_splice_site(pos_in_window: int, splice_track: Optional[np.ndarray],
                       resolution: int = 128) -> float:
    """
    Return predicted splice-site score at a position from an AlphaGenome
    SPLICE_SITES track. Returns 0.5 (neutral) if no track is provided.
    """
    if splice_track is None:
        return 0.5
    bin_idx = pos_in_window // resolution
    if 0 <= bin_idx < len(splice_track):
        return float(splice_track[bin_idx])
    return 0.5


def predict_isoforms(
    genes_df: pd.DataFrame,
    ref_seq: str,
    mut_seq: str,
    win_start: int,
    edit_start: int,
    edit_end: int,
    mut_shift: int,
    ref_splice_track: Optional[np.ndarray] = None,
    mut_splice_track: Optional[np.ndarray] = None,
    splice_threshold: float = 0.2,
) -> Dict[str, Dict]:
    """
    Predict the mature mRNA(s) following DNA edits.

    When AlphaGenome SPLICE_SITES tracks are provided, exons whose donor or
    acceptor splice sites drop below `splice_threshold` in the mutant are
    flagged as potentially skipped.

    Args:
        genes_df: GTF annotations DataFrame from ReferenceContextBuilder.
        ref_seq / mut_seq: 1 Mb reference and mutant sequence windows.
        win_start: Genomic start of the window.
        edit_start / edit_end: Coordinates of the DNA modification.
        mut_shift: Net length change from the edit (negative for deletions).
        ref_splice_track: 1-D array of splice-site scores for reference (optional).
        mut_splice_track: 1-D array of splice-site scores for mutant (optional).
        splice_threshold: Mutant splice score below which an exon is flagged.

    Returns:
        Dict mapping transcript_id → {
            'ref_mrna': str,
            'mut_mrna': str,
            'skipped_exons': list of exon indices that may be skipped in mutant,
            'frameshift': bool,
            'splice_disrupted': bool,
        }
    """
    isoforms: Dict[str, Dict] = {}
    transcripts = get_transcript_ids(genes_df)

    for tid in transcripts:
        exons = extract_exons(genes_df, tid)
        if exons.empty:
            continue

        strand = exons.iloc[0]['Strand']
        ref_mrna_dna = ""
        mut_mrna_dna = ""
        skipped_exons: List[int] = []

        for exon_idx, (_, exon) in enumerate(exons.iterrows()):
            exon_start = int(exon['Start'])
            exon_end   = int(exon['End'])

            # ---- Reference mRNA assembly ----
            rel_start = max(exon_start - win_start, 0)
            rel_end   = min(exon_end   - win_start, len(ref_seq))
            if rel_start < rel_end:
                ref_mrna_dna += ref_seq[rel_start:rel_end]

            # ---- Check splice-site integrity in mutant ----
            # Check the acceptor site (exon start) and donor site (exon end)
            mut_acceptor_pos = get_mut_coord(exon_start, win_start, edit_start, edit_end, mut_shift)
            mut_donor_pos    = get_mut_coord(exon_end,   win_start, edit_start, edit_end, mut_shift)

            acceptor_score = _score_splice_site(mut_acceptor_pos, mut_splice_track)
            donor_score    = _score_splice_site(mut_donor_pos,    mut_splice_track)

            # If either splice site drops below threshold, flag the exon for skipping
            splice_disrupted_here = (
                mut_splice_track is not None
                and (acceptor_score < splice_threshold or donor_score < splice_threshold)
            )

            # ---- Mutant mRNA assembly ----
            mut_start = max(mut_acceptor_pos, 0)
            mut_end   = min(mut_donor_pos,    len(mut_seq))
            if mut_start < mut_end and not splice_disrupted_here:
                mut_mrna_dna += mut_seq[mut_start:mut_end]
            elif splice_disrupted_here:
                skipped_exons.append(exon_idx)

        # Reverse complement for minus-strand genes
        if strand == '-':
            ref_mrna_dna = str(Seq(ref_mrna_dna).reverse_complement())
            mut_mrna_dna = str(Seq(mut_mrna_dna).reverse_complement())

        ref_mrna = ref_mrna_dna.upper().replace('T', 'U')
        mut_mrna = mut_mrna_dna.upper().replace('T', 'U')

        # Frameshift detection: length difference not divisible by 3
        len_diff = len(mut_mrna) - len(ref_mrna)
        frameshift = (len_diff % 3 != 0) and (len(ref_mrna) > 0)

        isoforms[tid] = {
            'ref_mrna': ref_mrna,
            'mut_mrna': mut_mrna,
            'skipped_exons': skipped_exons,
            'frameshift': frameshift,
            'splice_disrupted': len(skipped_exons) > 0,
        }

    return isoforms

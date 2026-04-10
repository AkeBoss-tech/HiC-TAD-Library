from Bio.Seq import Seq
from typing import Any, Dict

def find_longest_orf(mrna_seq: str) -> str:
    """
    Finds the longest Open Reading Frame (ORF) in the given mRNA sequence.
    An ORF starts with 'AUG' and ends with a stop codon ('UAA', 'UAG', 'UGA').
    """
    longest_orf = ""
    start_codon = "AUG"
    stop_codons = {"UAA", "UAG", "UGA"}

    for i in range(len(mrna_seq) - 2):
        if mrna_seq[i:i+3] == start_codon:
            for j in range(i + 3, len(mrna_seq) - 2, 3):
                codon = mrna_seq[j:j+3]
                if codon in stop_codons:
                    orf = mrna_seq[i:j+3]
                    if len(orf) > len(longest_orf):
                        longest_orf = orf
                    break

    return longest_orf

def translate(mrna_seq: str) -> str:
    """
    Translates an mRNA sequence into an amino acid sequence.
    Extracts the longest ORF and translates it.
    """
    orf = find_longest_orf(mrna_seq)
    if not orf:
        return ""

    # Biopython's translate function works with standard RNA alphabets
    seq_obj = Seq(orf)
    # to_stop=True ensures we stop translating at the first stop codon
    return str(seq_obj.translate(to_stop=True))

def compare_translation(ref_mrna: str, mut_mrna: str) -> Dict[str, Any]:
    """
    Compares the translation of reference and mutant mRNA sequences.
    Detects frameshifts or truncation events.
    """
    ref_aa = translate(ref_mrna)
    mut_aa = translate(mut_mrna)

    frameshift = False

    if len(ref_aa) != len(mut_aa) and len(ref_aa) > 0:
        # A simple heuristic for a frameshift: the sequence lengths differ
        # and the sequences aren't identical up to the length of the shorter one
        min_len = min(len(ref_aa), len(mut_aa))
        if ref_aa[:min_len] != mut_aa[:min_len]:
            frameshift = True
        elif len(mut_aa) < len(ref_aa):
            # It could also be a frameshift leading to premature stop codon
            # Or a simple truncation
            frameshift = True

    return {
        'ref_aa': ref_aa,
        'mut_aa': mut_aa,
        'is_frameshift': frameshift
    }

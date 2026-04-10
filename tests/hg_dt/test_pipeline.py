import pytest
import pandas as pd
from src.hg_dt.translate.transcript import predict_isoforms
from src.hg_dt.translate.translator import compare_translation, translate
from src.hg_dt.models.protein import ProteinFolder
from src.hg_dt.analyze.deltas import compute_structure_impact

def test_frameshift_deletion():
    # A simple gene setup with 1 exon.
    # Ref sequence window
    ref_seq = "ATGGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    # 1 bp deletion at position 4 (G -> deleted)
    # Ref: ATG GCC ATT GTA ATG GGC CGC TGA (start -> GCC -> ATT -> GTA -> ATG -> GGC -> CGC -> TGA stop)
    # Ref AA: M A I V M G R *

    # Mut: ATG CCA TTG TAA (start -> CCA -> TTG -> TAA stop)
    # Mut AA: M P L *

    mut_seq = "ATGCCATTGTAATGGGCCGCTGAAAGGGTGCCCGATAG"

    genes_data = {
        'Feature': ['exon'],
        'transcript_id': ['T1'],
        'Start': [0],
        'End': [len(ref_seq)],
        'Strand': ['+']
    }
    genes_df = pd.DataFrame(genes_data)

    # Predict isoforms
    # edit_start = 3, edit_end = 4, mut_shift = -1
    isoforms = predict_isoforms(genes_df, ref_seq, mut_seq, win_start=0, edit_start=3, edit_end=4, mut_shift=-1)

    assert 'T1' in isoforms
    t1_iso = isoforms['T1']

    ref_mrna = t1_iso['ref_mrna']
    mut_mrna = t1_iso['mut_mrna']

    # Check translation
    comp = compare_translation(ref_mrna, mut_mrna)
    assert comp['is_frameshift'] is True

    # Check AA sequences
    assert comp['ref_aa'] == "MAIVMGR"
    assert comp['mut_aa'] == "MPL"

def test_structural_collapse():
    # Deterministic: avoid live ESMFold / mock randomness from ProteinFolder.
    ref_plddt = [92.0] * 60
    mut_plddt = [42.0] * 15
    impact = compute_structure_impact(ref_plddt, mut_plddt)
    assert impact['collapse'] is True
    assert impact['avg_drop'] > 0
    assert impact['avg_ref_plddt'] > 70.0
    assert impact['avg_mut_plddt'] < 50.0

def test_no_structural_collapse():
    folder = ProteinFolder()

    # WT
    ref_aa = "M" * 60
    # Mut with synonymous or non-truncating that still folds well
    mut_aa = "M" * 60

    ref_struct = folder.predict_structure(ref_aa)
    mut_struct = folder.predict_structure(mut_aa)

    impact = compute_structure_impact(ref_struct['plddt'], mut_struct['plddt'])

    assert impact['collapse'] is False
    assert impact['avg_drop'] == 0.0

import numpy as np

import pipeline.stage0_variant_context as stage0_variant_context
from pipeline.stage1_genomics import derive_variant_sequences
from pipeline.stage2_protein import cosine_distance
from pipeline.stage3_folding import kabsch_rmsd


def test_build_context_has_boundary_collapse_metadata(monkeypatch):
    monkeypatch.setattr(stage0_variant_context, "fetch_ensembl_isoforms", lambda: ([], None))
    context = stage0_variant_context.build_context()
    assert context["boundary_collapse"]["variant_id"] == "del_4ctcf"
    assert context["boundary_collapse"]["delta_insulation"] > 0.2
    assert "fallback_variant" in context


def test_derive_variant_sequences_uses_isoform_tail_before_proxy():
    wt_domain = "A" * 79
    full_seq = "M" * 900 + wt_domain
    context = {
        "isoforms": [
            {
                "transcript_id": "ENSMUST_TEST",
                "display_name": "Unc5b-test-isoform",
                "protein_sequence": "M" * 880 + ("B" * 79),
            }
        ],
        "fallback_variant": {"truncate_residues": 16},
    }

    records = derive_variant_sequences(full_seq, wt_domain, context)
    assert [record["id"] for record in records] == ["wt_canonical", "ENSMUST_TEST"]
    assert records[1]["sequence"] == "B" * 79


def test_derive_variant_sequences_falls_back_to_proxy():
    wt_domain = "A" * 79
    records = derive_variant_sequences("M" * 900 + wt_domain, wt_domain, {"isoforms": [], "fallback_variant": {"truncate_residues": 16}})
    assert records[1]["id"] == "boundary_collapse_proxy"
    assert len(records[1]["sequence"]) == 63


def test_cosine_distance_zero_for_identical_vectors():
    a = np.array([1.0, 2.0, 3.0])
    assert cosine_distance(a, a) == 0.0


def test_kabsch_rmsd_zero_for_rigid_transform():
    ref = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
    mob = ref + np.array([10.0, -3.0, 5.0])
    assert kabsch_rmsd(ref, mob) < 1e-6

import os
import numpy as np
import pytest

from src.hg_dt.analyze.deltas import (
    contact_delta,
    accessibility_delta,
    expression_delta,
    find_distal_loops,
    find_silenced_elements
)
from src.polymer_sim import polymer_from_contact_map
from src.tad_boundaries import compute_insulation_delta

# Mock AlphaGenomeConnector since we might not have a real API key in tests
class MockAlphaGenomeConnector:
    def predict_variant(self, chrom, start, end, variant_pos, ref_base, alt_base, organism="HUMAN", requested_outputs=None):
        class MockTrack:
            def __init__(self, data):
                self.data = data
                self.values = data

            def __array__(self):
                return self.values

        class MockOutput:
            def __init__(self, is_mutant=False):
                # Contact map (100x100)
                base_map = np.zeros((100, 100))
                # Add diagonal
                for i in range(100):
                    base_map[i, i] = 10.0
                    if i < 99:
                        base_map[i, i+1] = base_map[i+1, i] = 5.0

                # Add loop (e.g. enhancer to promoter at distance 40)
                if not is_mutant:
                    base_map[20, 60] = base_map[60, 20] = 2.0  # MYC enhancer loop

                self.contact_maps = MockTrack(base_map)

                # ATAC track
                atac = np.ones(100) * 0.1
                atac[20] = 5.0 # enhancer
                atac[60] = 3.0 # promoter
                if is_mutant:
                    atac[20] = 0.5 # enhancer silenced

                self.atac = MockTrack(atac)

                # Expression track (CAGE)
                expr = np.zeros(100)
                expr[60] = 10.0 # expression at promoter
                if is_mutant:
                    expr[60] = 2.0 # expression dropped

                self.cage = MockTrack(expr)
                self.rna_seq = MockTrack(expr)

            def get(self, track_type):
                if track_type.name == 'CONTACT_MAPS':
                    return self.contact_maps
                elif track_type.name == 'ATAC':
                    return self.atac
                elif track_type.name == 'CAGE':
                    return self.cage
                return None

        class MockVariantOutput:
            def __init__(self):
                self.reference = MockOutput(is_mutant=False)
                self.alternate = MockOutput(is_mutant=True)

        return MockVariantOutput()


def test_myc_enhancer_deletion():
    """Integration test for MYC enhancer deletion."""

    # Use mock client
    client = MockAlphaGenomeConnector()

    outputs = client.predict_variant(
        chrom='chr8',
        start=127000000,
        end=128000000,
        variant_pos=127500000, # mock enhancer
        ref_base='A',
        alt_base='T'
    )

    # 1. Extract tracks
    ref_map = np.array(outputs.reference.contact_maps)
    mut_map = np.array(outputs.alternate.contact_maps)

    ref_atac = np.array(outputs.reference.atac)
    mut_atac = np.array(outputs.alternate.atac)

    ref_expr = np.array(outputs.reference.cage)
    mut_expr = np.array(outputs.alternate.cage)

    # 2. Compute deltas
    c_delta = contact_delta(ref_map, mut_map)
    a_delta = accessibility_delta(ref_atac, mut_atac)
    e_delta = expression_delta(ref_expr, mut_expr)

    # 3. Find loops in reference
    ref_loops = find_distal_loops(ref_map, min_dist=30, min_freq=1.0)

    # The enhancer loop (20, 60) should be found in ref
    assert (20, 60) in ref_loops

    # The loop should be absent or weakened in mutant
    mut_loops = find_distal_loops(mut_map, min_dist=30, min_freq=1.0)
    assert (20, 60) not in mut_loops

    # 4. Verify contact delta shows loss of loop
    assert c_delta[20, 60] < -1.0 # significant loss

    # 5. Verify silenced elements
    silenced = find_silenced_elements(ref_atac, mut_atac, threshold=0.5)
    assert 20 in silenced # enhancer at index 20 silenced

    # 6. Verify "Likelihood of Expression Loss"
    # For target gene at index 60
    expr_loss_log2 = e_delta[60]
    assert expr_loss_log2 < -1.0 # more than 2x drop

    likelihood = 1.0 - (2 ** expr_loss_log2)
    assert likelihood > 0.5 # greater than 50% likelihood of loss

    # 7. Check polymer sim runs
    coords_ref = polymer_from_contact_map(ref_map, n_steps=100)
    coords_mut = polymer_from_contact_map(mut_map, n_steps=100)

    assert coords_ref.shape == (100, 3)
    assert coords_mut.shape == (100, 3)

"""
Tests for src/polymer_sim.py - 3D polymer simulation from contact matrices.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from src.polymer_sim import (
    contact_matrix_to_restraints,
    insulation_to_backbone_stiffness,
    simulate_polymer,
    polymer_from_cooler,
    _init_coords,
)


class TestContactMatrixToRestraints:
    """Test suite for contact matrix to restraints conversion."""

    @pytest.mark.unit
    def test_basic_restraints_generation(self, mock_contact_matrix):
        """Test that restraints are generated from contact matrix."""
        restraints = contact_matrix_to_restraints(
            mock_contact_matrix,
            contact_threshold_quantile=0.70
        )

        assert isinstance(restraints, list)
        assert len(restraints) > 0

        # Check restraint format: (i, j, rest_length, k)
        for r in restraints:
            assert len(r) == 4
            i, j, rest_len, k = r
            assert isinstance(i, int)
            assert isinstance(j, int)
            assert isinstance(rest_len, float)
            assert isinstance(k, float)
            assert i < j  # Should only have upper triangle

    @pytest.mark.unit
    def test_restraints_respect_threshold(self, mock_contact_matrix):
        """Test that threshold affects number of restraints."""
        low_threshold = contact_matrix_to_restraints(
            mock_contact_matrix,
            contact_threshold_quantile=0.5
        )

        high_threshold = contact_matrix_to_restraints(
            mock_contact_matrix,
            contact_threshold_quantile=0.9
        )

        # Lower threshold should produce more restraints
        assert len(low_threshold) > len(high_threshold)

    @pytest.mark.unit
    def test_restraints_spring_constants_in_range(self, mock_contact_matrix):
        """Test that spring constants are within specified range."""
        k_min, k_max = 0.5, 8.0
        restraints = contact_matrix_to_restraints(
            mock_contact_matrix,
            contact_threshold_quantile=0.70,
            k_min=k_min,
            k_max=k_max
        )

        for i, j, rest_len, k in restraints:
            assert k_min <= k <= k_max

    @pytest.mark.unit
    def test_restraints_rest_length_range(self, mock_contact_matrix):
        """Test that rest lengths are in expected range."""
        restraints = contact_matrix_to_restraints(
            mock_contact_matrix,
            contact_threshold_quantile=0.70
        )

        for i, j, rest_len, k in restraints:
            # Default range: 1.0 - 3.0
            assert 1.0 <= rest_len <= 3.0

    @pytest.mark.unit
    def test_empty_matrix_handling(self):
        """Test handling of empty or all-zero matrix."""
        empty_matrix = np.zeros((50, 50))
        restraints = contact_matrix_to_restraints(empty_matrix)

        assert isinstance(restraints, list)
        assert len(restraints) == 0

    @pytest.mark.unit
    def test_nan_handling(self):
        """Test that NaN values in matrix are handled correctly."""
        matrix = np.random.random((50, 50))
        matrix[10:20, 10:20] = np.nan

        restraints = contact_matrix_to_restraints(matrix)

        # Should still produce valid restraints
        assert isinstance(restraints, list)
        for i, j, rest_len, k in restraints:
            assert not np.isnan(rest_len)
            assert not np.isnan(k)

    @pytest.mark.unit
    def test_symmetric_matrix_handling(self):
        """Test that asymmetric matrices are symmetrized."""
        n = 50
        matrix = np.random.random((n, n))

        restraints = contact_matrix_to_restraints(matrix)

        # Should only have upper triangle restraints (no duplicates)
        pairs = set()
        for i, j, _, _ in restraints:
            pair = (min(i, j), max(i, j))
            assert pair not in pairs
            pairs.add(pair)


class TestInsulationToBackbone:
    """Test suite for insulation to backbone stiffness mapping."""

    @pytest.mark.unit
    def test_backbone_stiffness_output_length(self, mock_insulation_scores):
        """Test that output has correct length (N-1 for N bins)."""
        n_bins = len(mock_insulation_scores)
        stiffness = insulation_to_backbone_stiffness(mock_insulation_scores)

        assert len(stiffness) == n_bins - 1

    @pytest.mark.unit
    def test_backbone_stiffness_range(self, mock_insulation_scores):
        """Test that stiffness values are within specified range."""
        k_soft, k_stiff = 2.0, 15.0
        stiffness = insulation_to_backbone_stiffness(
            mock_insulation_scores,
            k_soft=k_soft,
            k_stiff=k_stiff
        )

        assert (stiffness >= k_soft).all()
        assert (stiffness <= k_stiff).all()

    @pytest.mark.unit
    def test_low_insulation_high_stiffness(self):
        """Test that low insulation (boundaries) produces high stiffness."""
        # Create scores with clear boundary
        scores = np.zeros(50)
        scores[24:26] = -2.0  # Strong boundary

        stiffness = insulation_to_backbone_stiffness(
            scores,
            k_soft=1.0,
            k_stiff=10.0
        )

        # Stiffness at boundary should be high
        boundary_idx = 24
        # Check neighborhood around boundary
        assert stiffness[boundary_idx] > np.mean(stiffness)

    @pytest.mark.unit
    def test_uniform_insulation(self):
        """Test handling of uniform insulation scores."""
        uniform_scores = np.ones(50) * 0.5
        stiffness = insulation_to_backbone_stiffness(uniform_scores)

        # Should return mid-range values
        assert np.allclose(stiffness, stiffness[0])

    @pytest.mark.unit
    def test_nan_handling_in_insulation(self):
        """Test that NaN values in insulation are handled."""
        scores = np.random.randn(50)
        scores[10:15] = np.nan

        stiffness = insulation_to_backbone_stiffness(scores)

        assert len(stiffness) == 49
        assert not np.any(np.isnan(stiffness))


class TestSimulatePolymer:
    """Test suite for polymer simulation."""

    @pytest.mark.unit
    def test_simulate_polymer_output_shape(self, sample_restraints):
        """Test that simulation output has correct shape."""
        n_beads = 50
        coords = simulate_polymer(
            n_beads,
            sample_restraints,
            n_steps=100,
            seed=42
        )

        assert coords.shape == (n_beads, 3)

    @pytest.mark.unit
    def test_simulate_polymer_reproducibility(self, sample_restraints):
        """Test that simulation is reproducible with same seed."""
        n_beads = 30

        coords1 = simulate_polymer(n_beads, sample_restraints, n_steps=100, seed=42)
        coords2 = simulate_polymer(n_beads, sample_restraints, n_steps=100, seed=42)

        assert np.allclose(coords1, coords2)

    @pytest.mark.unit
    def test_simulate_polymer_different_seeds(self, sample_restraints):
        """Test that different seeds produce different results."""
        n_beads = 30

        coords1 = simulate_polymer(n_beads, sample_restraints, n_steps=100, seed=42)
        coords2 = simulate_polymer(n_beads, sample_restraints, n_steps=100, seed=123)

        assert not np.allclose(coords1, coords2)

    @pytest.mark.unit
    def test_simulate_polymer_no_restraints(self):
        """Test simulation with no Hi-C restraints (backbone only)."""
        n_beads = 30
        coords = simulate_polymer(
            n_beads,
            restraints=[],
            n_steps=100,
            seed=42
        )

        assert coords.shape == (n_beads, 3)
        # Check that beads are not all at the same position
        assert not np.allclose(coords, coords[0])

    @pytest.mark.unit
    def test_simulate_polymer_with_backbone_stiffness(self, sample_restraints):
        """Test simulation with variable backbone stiffness."""
        n_beads = 30
        backbone_k = np.linspace(1.0, 10.0, n_beads - 1)

        coords = simulate_polymer(
            n_beads,
            sample_restraints,
            backbone_k=backbone_k,
            n_steps=100,
            seed=42
        )

        assert coords.shape == (n_beads, 3)

    @pytest.mark.unit
    def test_simulate_polymer_chain_connectivity(self, sample_restraints):
        """Test that polymer chain remains connected."""
        n_beads = 30
        coords = simulate_polymer(
            n_beads,
            sample_restraints,
            n_steps=500,
            seed=42
        )

        # Check that consecutive beads are reasonably close
        distances = np.linalg.norm(coords[1:] - coords[:-1], axis=1)

        # All consecutive beads should be within reasonable distance
        # (backbone rest length + some tolerance)
        assert (distances < 5.0).all()

    @pytest.mark.unit
    def test_init_coords_shape(self):
        """Test coordinate initialization."""
        n_beads = 50
        rng = np.random.default_rng(42)
        coords = _init_coords(n_beads, rng)

        assert coords.shape == (n_beads, 3)
        # First bead should be at origin
        assert np.allclose(coords[0], [0, 0, 0])

    @pytest.mark.unit
    @pytest.mark.slow
    def test_simulate_polymer_convergence(self, sample_restraints):
        """Test that simulation converges with enough steps."""
        n_beads = 30

        # Short simulation
        coords_short = simulate_polymer(
            n_beads, sample_restraints, n_steps=100, seed=42
        )

        # Long simulation
        coords_long = simulate_polymer(
            n_beads, sample_restraints, n_steps=5000, seed=42
        )

        # Both should produce valid structures (not checking convergence to same point)
        assert coords_short.shape == coords_long.shape == (n_beads, 3)

    @pytest.mark.unit
    def test_simulate_polymer_temperature_effect(self, sample_restraints):
        """Test that temperature affects simulation."""
        n_beads = 30

        coords_cold = simulate_polymer(
            n_beads, sample_restraints, n_steps=500, temperature=0.1, seed=42
        )

        coords_hot = simulate_polymer(
            n_beads, sample_restraints, n_steps=500, temperature=5.0, seed=42
        )

        # Different temperatures should produce different structures
        assert not np.allclose(coords_cold, coords_hot)


class TestPolymerFromCooler:
    """Test suite for end-to-end cooler to polymer conversion."""

    @pytest.mark.unit
    def test_polymer_from_cooler_output_shape(self, mock_cooler):
        """Test that polymer_from_cooler produces correct output."""
        coords = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            n_steps=100,
            seed=42
        )

        expected_bins = (26_100_000 - 26_000_000) // 5000
        assert coords.shape == (expected_bins, 3)

    @pytest.mark.unit
    def test_polymer_from_cooler_with_insulation(self, mock_cooler, mock_insulation_scores):
        """Test polymer simulation with insulation scores."""
        # Need to match the region size
        region_bins = 20
        insulation = mock_insulation_scores[:region_bins]

        coords = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            insulation_scores=insulation,
            n_steps=100,
            seed=42
        )

        assert coords.shape[0] == region_bins
        assert coords.shape[1] == 3

    @pytest.mark.unit
    def test_polymer_from_cooler_reproducibility(self, mock_cooler):
        """Test reproducibility of polymer_from_cooler."""
        coords1 = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            n_steps=100,
            seed=42
        )

        coords2 = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            n_steps=100,
            seed=42
        )

        assert np.allclose(coords1, coords2)

    @pytest.mark.unit
    def test_polymer_from_cooler_threshold_parameter(self, mock_cooler):
        """Test that contact threshold parameter is used."""
        coords_low = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            contact_threshold_quantile=0.5,
            n_steps=100,
            seed=42
        )

        coords_high = polymer_from_cooler(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            contact_threshold_quantile=0.9,
            n_steps=100,
            seed=42
        )

        # Both should produce valid structures
        assert coords_low.shape == coords_high.shape
        # Different thresholds should affect the structure
        assert not np.allclose(coords_low, coords_high)


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.mark.unit
    def test_small_matrix(self):
        """Test with very small contact matrix."""
        small_matrix = np.array([
            [1.0, 0.5, 0.2],
            [0.5, 1.0, 0.5],
            [0.2, 0.5, 1.0]
        ])

        restraints = contact_matrix_to_restraints(small_matrix)
        coords = simulate_polymer(3, restraints, n_steps=50)

        assert coords.shape == (3, 3)

    @pytest.mark.unit
    def test_large_n_beads(self):
        """Test with larger number of beads."""
        n_beads = 500
        restraints = [(i, i+10, 2.0, 1.0) for i in range(0, n_beads-10, 20)]

        coords = simulate_polymer(
            n_beads,
            restraints,
            n_steps=100,
            seed=42
        )

        assert coords.shape == (n_beads, 3)

    @pytest.mark.unit
    def test_invalid_restraints_handling(self):
        """Test that invalid restraints are handled gracefully."""
        n_beads = 30
        # Restraints with invalid indices should not crash
        invalid_restraints = [
            (0, 5, 2.0, 1.0),
            (100, 105, 2.0, 1.0),  # Out of bounds
        ]

        # This might raise an error or skip invalid restraints
        # depending on implementation
        try:
            coords = simulate_polymer(n_beads, invalid_restraints, n_steps=10)
            # If it succeeds, check output is valid
            assert coords.shape == (n_beads, 3)
        except (IndexError, ValueError):
            # It's also acceptable to raise an error
            pass

    @pytest.mark.unit
    def test_zero_steps(self):
        """Test simulation with zero steps returns initialized coords."""
        n_beads = 20
        coords = simulate_polymer(
            n_beads,
            restraints=[],
            n_steps=0,
            seed=42
        )

        assert coords.shape == (n_beads, 3)

    @pytest.mark.unit
    def test_negative_spring_constant_handling(self):
        """Test handling of edge case with k_min/k_max."""
        matrix = np.random.random((30, 30))

        # Should clamp to valid range
        restraints = contact_matrix_to_restraints(
            matrix,
            k_min=0.0,
            k_max=10.0
        )

        for i, j, rest_len, k in restraints:
            assert k >= 0.0
            assert k <= 10.0


class TestNumericalStability:
    """Test numerical stability of simulation."""

    @pytest.mark.unit
    def test_no_nan_in_output(self, sample_restraints):
        """Test that simulation doesn't produce NaN values."""
        coords = simulate_polymer(
            50,
            sample_restraints,
            n_steps=500,
            seed=42
        )

        assert not np.any(np.isnan(coords))
        assert not np.any(np.isinf(coords))

    @pytest.mark.unit
    def test_no_explosion(self, sample_restraints):
        """Test that coordinates don't explode to large values."""
        coords = simulate_polymer(
            50,
            sample_restraints,
            n_steps=1000,
            seed=42
        )

        # Coordinates should be within reasonable range
        # (not checking specific values, just no explosion)
        assert np.abs(coords).max() < 1000.0

    @pytest.mark.unit
    def test_extreme_temperature(self, sample_restraints):
        """Test simulation with extreme temperature values."""
        # Very low temperature
        coords_cold = simulate_polymer(
            30, sample_restraints, n_steps=100, temperature=0.01, seed=42
        )
        assert not np.any(np.isnan(coords_cold))

        # High temperature (but not too extreme)
        coords_hot = simulate_polymer(
            30, sample_restraints, n_steps=100, temperature=10.0, seed=42
        )
        assert not np.any(np.isnan(coords_hot))

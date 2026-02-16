"""
Tests for src/tad_boundaries.py - TAD boundary detection and analysis.
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch

from src.tad_boundaries import (
    parse_coordinates,
    make_view_df,
    compute_directionality_index,
    score_boundary_prominence,
    call_tad_intervals,
    multiscale_insulation,
    classify_boundary_persistence,
    boundary_pileup,
    _li_threshold,
)


class TestHelperFunctions:
    """Test suite for helper functions."""

    @pytest.mark.unit
    def test_parse_coordinates_basic(self):
        """Test parsing of basic coordinate strings."""
        chrom, start, end = parse_coordinates("chr12:26,000,000-28,000,000")
        assert chrom == "chr12"
        assert start == 26_000_000
        assert end == 28_000_000

    @pytest.mark.unit
    def test_parse_coordinates_no_commas(self):
        """Test parsing coordinates without comma separators."""
        chrom, start, end = parse_coordinates("chr13:83500000-84500000")
        assert chrom == "chr13"
        assert start == 83_500_000
        assert end == 84_500_000

    @pytest.mark.unit
    def test_parse_coordinates_different_chromosomes(self):
        """Test parsing with various chromosome names."""
        test_cases = [
            ("chr1:1,000-2,000", "chr1", 1000, 2000),
            ("chr2:100000-200000", "chr2", 100000, 200000),
            ("chrX:1,000,000-2,000,000", "chrX", 1_000_000, 2_000_000),
        ]
        for coord_str, exp_chr, exp_start, exp_end in test_cases:
            chrom, start, end = parse_coordinates(coord_str)
            assert chrom == exp_chr
            assert start == exp_start
            assert end == exp_end

    @pytest.mark.unit
    def test_make_view_df_structure(self):
        """Test that make_view_df creates correct DataFrame structure."""
        df = make_view_df("chr12", 26_000_000, 28_000_000)

        assert isinstance(df, pd.DataFrame)
        assert len(df) == 1
        assert list(df.columns) == ['chrom', 'start', 'end', 'name']
        assert df.iloc[0]['chrom'] == 'chr12'
        assert df.iloc[0]['start'] == 26_000_000
        assert df.iloc[0]['end'] == 28_000_000
        assert df.iloc[0]['name'] == 'chr12'

    @pytest.mark.unit
    def test_make_view_df_different_chromosomes(self):
        """Test make_view_df with different chromosome inputs."""
        df1 = make_view_df("chr1", 0, 1_000_000)
        df2 = make_view_df("chrX", 5_000_000, 10_000_000)

        assert df1.iloc[0]['chrom'] == 'chr1'
        assert df1.iloc[0]['name'] == 'chr1'
        assert df2.iloc[0]['chrom'] == 'chrX'


class TestDirectionalityIndex:
    """Test suite for Directionality Index computation."""

    @pytest.mark.unit
    def test_compute_directionality_index_output_structure(self, mock_cooler):
        """Test that DI output has correct structure."""
        result = compute_directionality_index(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            window_bp=50_000
        )

        assert isinstance(result, pd.DataFrame)
        assert 'chrom' in result.columns
        assert 'start' in result.columns
        assert 'end' in result.columns
        assert 'upstream_sum' in result.columns
        assert 'downstream_sum' in result.columns
        assert 'DI' in result.columns

    @pytest.mark.unit
    def test_compute_directionality_index_length(self, mock_cooler):
        """Test that DI output length matches region."""
        result = compute_directionality_index(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            window_bp=25_000
        )

        expected_bins = (26_100_000 - 26_000_000) // 5000
        assert len(result) == expected_bins

    @pytest.mark.unit
    def test_directionality_index_values_are_numeric(self, mock_cooler):
        """Test that DI values are numeric (may contain NaN at edges)."""
        result = compute_directionality_index(
            mock_cooler,
            "chr12:26,000,000-26,100,000",
            window_bp=25_000
        )

        assert result['DI'].dtype in [np.float64, np.float32, float]
        assert result['upstream_sum'].dtype in [np.float64, np.float32, float]
        assert result['downstream_sum'].dtype in [np.float64, np.float32, float]


class TestBoundaryProminence:
    """Test suite for boundary prominence scoring."""

    @pytest.mark.unit
    def test_score_boundary_prominence_output_structure(self, sample_insulation_table):
        """Test that prominence scoring output has correct structure."""
        result = score_boundary_prominence(
            sample_insulation_table,
            window_bp=25000,
            prominence_thresholds=(0.2, 0.5)
        )

        assert isinstance(result, pd.DataFrame)
        assert 'prominence' in result.columns
        assert 'boundary_class' in result.columns
        assert len(result) == len(sample_insulation_table)

    @pytest.mark.unit
    def test_score_boundary_prominence_classes(self, sample_insulation_table):
        """Test that boundary classes are assigned correctly."""
        result = score_boundary_prominence(
            sample_insulation_table,
            window_bp=25000,
            prominence_thresholds=(0.2, 0.5)
        )

        # Check that valid classes are assigned
        valid_classes = {'strong', 'weak', 'sub_threshold', None}
        unique_classes = set(result['boundary_class'].unique())
        assert unique_classes.issubset(valid_classes)

    @pytest.mark.unit
    def test_score_boundary_prominence_thresholds(self, sample_insulation_table):
        """Test that prominence thresholds work correctly."""
        result = score_boundary_prominence(
            sample_insulation_table,
            window_bp=25000,
            prominence_thresholds=(0.2, 0.5),
            min_distance_bins=3
        )

        # Boundaries should be detected
        boundaries = result[result['boundary_class'].notna()]
        assert len(boundaries) > 0

        # Check that prominence values match class assignments
        strong = result[result['boundary_class'] == 'strong']
        if len(strong) > 0:
            assert (strong['prominence'] >= 0.5).all()

        weak = result[result['boundary_class'] == 'weak']
        if len(weak) > 0:
            assert (weak['prominence'] >= 0.2).all()
            assert (weak['prominence'] < 0.5).all()

    @pytest.mark.unit
    def test_li_threshold_basic(self):
        """Test Li's thresholding algorithm."""
        # Create bimodal distribution
        values = np.concatenate([
            np.random.normal(2.0, 0.5, 100),
            np.random.normal(6.0, 0.5, 100)
        ])

        threshold = _li_threshold(values)

        # Threshold should be between the two modes
        assert 3.0 < threshold < 5.0

    @pytest.mark.unit
    def test_li_threshold_empty_array(self):
        """Test Li's threshold with empty array."""
        values = np.array([])
        threshold = _li_threshold(values)
        assert threshold == 0.0

    @pytest.mark.unit
    def test_li_threshold_with_nans(self):
        """Test Li's threshold handles NaN values."""
        values = np.array([1.0, 2.0, np.nan, 3.0, np.nan, 4.0])
        threshold = _li_threshold(values)
        assert not np.isnan(threshold)


class TestTADIntervals:
    """Test suite for TAD interval calling."""

    @pytest.mark.unit
    def test_call_tad_intervals_output_structure(self, sample_boundary_df, mock_cooler):
        """Test that TAD intervals have correct structure."""
        result = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='weak'
        )

        assert isinstance(result, pd.DataFrame)
        expected_cols = [
            'chrom', 'start', 'end', 'tad_id',
            'length_bp', 'n_bins',
            'left_boundary_class', 'right_boundary_class'
        ]
        for col in expected_cols:
            assert col in result.columns

    @pytest.mark.unit
    def test_call_tad_intervals_creates_intervals(self, sample_boundary_df, mock_cooler):
        """Test that TAD intervals are created from boundaries."""
        result = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='weak'
        )

        # Should have at least some TADs
        assert len(result) > 0

        # Each TAD should have positive length
        assert (result['length_bp'] > 0).all()
        assert (result['n_bins'] > 0).all()

    @pytest.mark.unit
    def test_call_tad_intervals_respects_class_filter(self, sample_boundary_df, mock_cooler):
        """Test that boundary class filtering works."""
        # Get TADs using only strong boundaries
        strong_only = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='strong'
        )

        # Get TADs using weak or strong
        weak_and_strong = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='weak'
        )

        # Should have same or fewer TADs with strong only
        assert len(strong_only) <= len(weak_and_strong)

    @pytest.mark.unit
    def test_call_tad_intervals_max_length_filter(self, sample_boundary_df, mock_cooler):
        """Test that maximum TAD length filter works."""
        result = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='weak',
            max_tad_length_bp=100_000
        )

        # All TADs should be below max length
        if len(result) > 0:
            assert (result['length_bp'] <= 100_000).all()

    @pytest.mark.unit
    def test_call_tad_intervals_tad_ids_unique(self, sample_boundary_df, mock_cooler):
        """Test that TAD IDs are unique."""
        result = call_tad_intervals(
            sample_boundary_df,
            mock_cooler,
            boundary_class_min='weak'
        )

        if len(result) > 0:
            assert len(result['tad_id'].unique()) == len(result)


class TestMultiscaleInsulation:
    """Test suite for multi-scale insulation analysis."""

    @pytest.mark.unit
    def test_multiscale_insulation_output_structure(self, mock_cooler):
        """Test multi-scale insulation output structure."""
        with patch('src.tad_boundaries.cooltools.insulation') as mock_ins:
            # Create mock insulation output
            mock_df = pd.DataFrame({
                'chrom': ['chr12'] * 100,
                'start': range(26_000_000, 26_500_000, 5000),
                'end': range(26_005_000, 26_505_000, 5000),
                'log2_insulation_score_25000': np.random.randn(100),
                'log2_insulation_score_50000': np.random.randn(100),
            })
            mock_ins.return_value = mock_df

            result, windows = multiscale_insulation(
                mock_cooler,
                "chr12:26,000,000-26,500,000",
                window_sizes_bp=[25000, 50000]
            )

            assert isinstance(result, pd.DataFrame)
            assert isinstance(windows, list)
            assert len(result) > 0

    @pytest.mark.unit
    def test_multiscale_insulation_auto_windows(self, mock_cooler):
        """Test automatic window size generation."""
        with patch('src.tad_boundaries.cooltools.insulation') as mock_ins:
            mock_df = pd.DataFrame({
                'chrom': ['chr12'] * 100,
                'start': range(26_000_000, 26_500_000, 5000),
                'end': range(26_005_000, 26_505_000, 5000),
            })
            # Add columns for each window
            for w in [25000, 50000, 100000]:
                mock_df[f'log2_insulation_score_{w}'] = np.random.randn(100)

            mock_ins.return_value = mock_df

            result, windows = multiscale_insulation(
                mock_cooler,
                "chr12:26,000,000-26,500,000",
                window_sizes_bp=None,
                n_windows=5
            )

            # Should generate multiple windows
            assert len(windows) > 1
            # Windows should be sorted
            assert windows == sorted(windows)


class TestBoundaryPersistence:
    """Test suite for boundary persistence classification."""

    @pytest.mark.unit
    def test_classify_boundary_persistence_output(self, sample_insulation_table):
        """Test boundary persistence output structure."""
        # Add multiple window columns
        for w in [25000, 50000, 100000]:
            if f'log2_insulation_score_{w}' not in sample_insulation_table.columns:
                sample_insulation_table[f'log2_insulation_score_{w}'] = np.random.randn(
                    len(sample_insulation_table)
                )

        result = classify_boundary_persistence(
            sample_insulation_table,
            window_sizes_bp=[25000, 50000, 100000],
            prominence_threshold=0.2
        )

        assert 'n_scales_boundary' in result.columns
        assert 'fraction_scales' in result.columns
        assert 'boundary_type' in result.columns
        assert 'min_scale_bp' in result.columns
        assert 'max_scale_bp' in result.columns

    @pytest.mark.unit
    def test_boundary_persistence_types(self, sample_insulation_table):
        """Test that valid boundary types are assigned."""
        for w in [25000, 50000, 100000]:
            sample_insulation_table[f'log2_insulation_score_{w}'] = np.random.randn(
                len(sample_insulation_table)
            )

        result = classify_boundary_persistence(
            sample_insulation_table,
            window_sizes_bp=[25000, 50000, 100000]
        )

        valid_types = {'constitutive', 'nested_sub', 'nested_meta', 'mixed', None}
        unique_types = set(result['boundary_type'].unique())
        assert unique_types.issubset(valid_types)


class TestBoundaryPileup:
    """Test suite for boundary pileup analysis."""

    @pytest.mark.unit
    def test_boundary_pileup_output_shape(self, mock_cooler, sample_boundary_df):
        """Test that pileup output has correct shape."""
        flank_bins = 20

        with patch('src.tad_boundaries.cooltools.expected_cis') as mock_exp:
            # Mock expected values
            mock_exp_df = pd.DataFrame({
                'balanced.avg': np.ones(100) * 0.5
            })
            mock_exp.return_value = mock_exp_df

            pileup, n_boundaries = boundary_pileup(
                mock_cooler,
                sample_boundary_df,
                flank_bins=flank_bins,
                boundary_class_min='weak'
            )

            expected_size = 2 * flank_bins + 1
            assert pileup.shape == (expected_size, expected_size)
            assert isinstance(n_boundaries, int)
            assert n_boundaries >= 0

    @pytest.mark.unit
    def test_boundary_pileup_respects_class_filter(self, mock_cooler, sample_boundary_df):
        """Test that class filtering affects pileup."""
        with patch('src.tad_boundaries.cooltools.expected_cis') as mock_exp:
            mock_exp_df = pd.DataFrame({'balanced.avg': np.ones(100) * 0.5})
            mock_exp.return_value = mock_exp_df

            _, n_strong = boundary_pileup(
                mock_cooler,
                sample_boundary_df,
                flank_bins=10,
                boundary_class_min='strong'
            )

            _, n_weak = boundary_pileup(
                mock_cooler,
                sample_boundary_df,
                flank_bins=10,
                boundary_class_min='weak'
            )

            # Should stack same or more boundaries with weak threshold
            assert n_weak >= n_strong

    @pytest.mark.unit
    def test_boundary_pileup_normalization_flag(self, mock_cooler, sample_boundary_df):
        """Test that normalization flag affects output."""
        with patch('src.tad_boundaries.cooltools.expected_cis') as mock_exp:
            mock_exp_df = pd.DataFrame({'balanced.avg': np.ones(100) * 0.5})
            mock_exp.return_value = mock_exp_df

            pileup_norm, _ = boundary_pileup(
                mock_cooler,
                sample_boundary_df,
                flank_bins=10,
                normalize_by_expected=True
            )

            pileup_raw, _ = boundary_pileup(
                mock_cooler,
                sample_boundary_df,
                flank_bins=10,
                normalize_by_expected=False
            )

            # Both should have same shape
            assert pileup_norm.shape == pileup_raw.shape


class TestEdgeCases:
    """Test edge cases and error handling."""

    @pytest.mark.unit
    def test_empty_boundary_df(self, mock_cooler):
        """Test handling of empty boundary DataFrame."""
        empty_df = pd.DataFrame({
            'chrom': [],
            'start': [],
            'end': [],
            'boundary_class': []
        })

        result = call_tad_intervals(empty_df, mock_cooler)

        # Should return empty DataFrame with correct columns
        assert len(result) == 0
        assert 'tad_id' in result.columns

    @pytest.mark.unit
    def test_single_boundary(self, mock_cooler):
        """Test handling of single boundary."""
        single_boundary = pd.DataFrame({
            'chrom': ['chr12'],
            'start': [26_200_000],
            'end': [26_205_000],
            'boundary_class': ['strong']
        })

        result = call_tad_intervals(single_boundary, mock_cooler)

        # Should create TADs on either side of the boundary
        assert len(result) >= 1

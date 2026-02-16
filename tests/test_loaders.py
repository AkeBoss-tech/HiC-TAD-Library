"""
Tests for src/loaders.py - Data loading utilities.
"""

import pytest
import os
from unittest.mock import patch, MagicMock
import cooler

from src.loaders import get_cooler, load_cooler


class TestGetCooler:
    """Test suite for get_cooler function."""

    @pytest.mark.unit
    def test_get_cooler_constructs_path_correctly(self):
        """Test that get_cooler constructs the correct file path."""
        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            get_cooler('test.mcool', resolution=5000)

            # Verify the URI was constructed correctly
            called_uri = mock_cooler_class.call_args[0][0]
            assert 'data/raw/test.mcool::resolutions/5000' in called_uri

    @pytest.mark.unit
    def test_get_cooler_mcool_uses_resolution(self):
        """Test that .mcool files use resolution in URI."""
        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            get_cooler('mouse_microc.mcool', resolution=10000)

            called_uri = mock_cooler_class.call_args[0][0]
            assert '::resolutions/10000' in called_uri

    @pytest.mark.unit
    def test_get_cooler_cool_no_resolution_suffix(self):
        """Test that .cool files don't add resolution suffix."""
        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            get_cooler('test.cool')

            called_uri = mock_cooler_class.call_args[0][0]
            assert '::resolutions' not in called_uri
            assert 'data/raw/test.cool' in called_uri

    @pytest.mark.unit
    def test_get_cooler_raises_on_missing_file(self):
        """Test that get_cooler raises FileNotFoundError for missing files."""
        with patch('src.loaders.os.path.exists', return_value=False):
            with pytest.raises(FileNotFoundError, match="File not found"):
                get_cooler('nonexistent.mcool')

    @pytest.mark.unit
    def test_get_cooler_default_resolution(self):
        """Test that default resolution is 1MB."""
        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            get_cooler('test.mcool')

            called_uri = mock_cooler_class.call_args[0][0]
            assert '::resolutions/1000000' in called_uri

    @pytest.mark.unit
    def test_load_cooler_alias(self):
        """Test that load_cooler is an alias for get_cooler."""
        assert load_cooler is get_cooler

    @pytest.mark.unit
    def test_get_cooler_absolute_path_construction(self):
        """Test that absolute paths are constructed correctly."""
        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.os.path.dirname') as mock_dirname, \
             patch('src.loaders.os.path.abspath') as mock_abspath, \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            # Mock the path construction
            mock_abspath.return_value = '/fake/path/src/loaders.py'
            mock_dirname.side_effect = ['/fake/path/src', '/fake/path']

            get_cooler('test.mcool', resolution=5000)

            # Verify that dirname was called correctly
            assert mock_dirname.call_count == 2

    @pytest.mark.unit
    def test_get_cooler_different_resolutions(self):
        """Test get_cooler with various common resolutions."""
        resolutions = [1000, 5000, 10000, 25000, 100000, 1000000]

        with patch('src.loaders.os.path.exists', return_value=True), \
             patch('src.loaders.cooler.Cooler') as mock_cooler_class:

            for res in resolutions:
                get_cooler('test.mcool', resolution=res)
                called_uri = mock_cooler_class.call_args[0][0]
                assert f'::resolutions/{res}' in called_uri


class TestLoaderIntegration:
    """Integration tests for loader functions."""

    @pytest.mark.integration
    @pytest.mark.requires_data
    def test_load_real_mcool_file_if_exists(self):
        """
        Test loading actual mcool file if it exists in data/raw.

        This test is skipped if the file doesn't exist.
        """
        # Construct expected path
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        test_file = os.path.join(base_dir, 'data', 'raw', 'mouse_microc.mcool')

        if not os.path.exists(test_file):
            pytest.skip("Test data file not found")

        # Attempt to load
        clr = get_cooler('mouse_microc.mcool', resolution=5000)

        # Verify it's a Cooler object
        assert isinstance(clr, cooler.Cooler)
        assert clr.binsize == 5000

        # Verify basic properties
        assert len(clr.chromsizes) > 0
        assert clr.info is not None

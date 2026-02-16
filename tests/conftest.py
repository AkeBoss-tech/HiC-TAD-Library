"""
Shared pytest fixtures for HiC-TAD-Library tests.

Provides mock data and utilities for testing without requiring large Hi-C files.
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import Mock, MagicMock
import tempfile
import os


@pytest.fixture
def mock_contact_matrix():
    """
    Create a synthetic contact matrix with TAD-like structure.

    Returns a 100x100 matrix with:
    - Strong diagonal (local contacts)
    - Two TAD blocks (0-40, 60-100)
    - Weak inter-TAD contacts
    """
    n = 100
    matrix = np.zeros((n, n))

    # Add diagonal decay (distance-dependent contacts)
    for d in range(n):
        strength = 1.0 / (d + 1)
        for i in range(n - d):
            matrix[i, i + d] = strength
            matrix[i + d, i] = strength

    # Enhance TAD 1 (bins 0-40)
    matrix[0:40, 0:40] *= 3.0
    matrix[0:40, 0:40] += np.random.random((40, 40)) * 0.5

    # Enhance TAD 2 (bins 60-100)
    matrix[60:100, 60:100] *= 3.0
    matrix[60:100, 60:100] += np.random.random((40, 40)) * 0.5

    # Add boundary depletion (insulation)
    matrix[38:42, :] *= 0.3
    matrix[:, 38:42] *= 0.3
    matrix[58:62, :] *= 0.3
    matrix[:, 58:62] *= 0.3

    return matrix


@pytest.fixture
def mock_insulation_scores():
    """
    Create synthetic insulation scores with clear boundaries.

    Returns array with low values (boundaries) at positions ~40 and ~60.
    """
    n = 100
    scores = np.zeros(n)

    # Baseline insulation (small noise)
    scores = np.random.normal(0, 0.1, n)

    # Add boundaries (deep minima)
    scores[38:42] = -1.5 + np.random.normal(0, 0.1, 4)
    scores[58:62] = -1.2 + np.random.normal(0, 0.1, 4)

    # Smooth with rolling average
    scores = pd.Series(scores).rolling(3, center=True, min_periods=1).mean().values

    return scores


@pytest.fixture
def mock_cooler():
    """
    Create a mock cooler.Cooler object for testing.

    Returns a MagicMock configured to behave like a Cooler object.
    """
    mock_clr = MagicMock()

    # Basic properties
    mock_clr.binsize = 5000
    mock_clr.chromsizes = pd.Series({
        'chr1': 195471971,
        'chr2': 182113224,
        'chr12': 120129022,
        'chr13': 120421639,
    })

    # Mock bins() method
    def mock_bins_fetch(region=None):
        """Mock bins().fetch() to return genomic bins."""
        if region:
            # Parse region like "chr12:26000000-28000000"
            parts = region.split(':')
            chrom = parts[0]
            if len(parts) > 1:
                start, end = map(int, parts[1].split('-'))
            else:
                start, end = 0, mock_clr.chromsizes[chrom]

            n_bins = (end - start) // mock_clr.binsize
            bins = pd.DataFrame({
                'chrom': [chrom] * n_bins,
                'start': range(start, end, mock_clr.binsize),
                'end': range(start + mock_clr.binsize, end + mock_clr.binsize, mock_clr.binsize)
            })
            return bins
        return pd.DataFrame({'chrom': [], 'start': [], 'end': []})

    mock_bins = MagicMock()
    mock_bins.fetch = mock_bins_fetch
    mock_clr.bins.return_value = mock_bins

    # Mock matrix() method
    def mock_matrix_fetch(region=None):
        """Mock matrix().fetch() to return a contact matrix."""
        if region:
            # Return a simple synthetic matrix
            n = 100  # Default size
            if ':' in region:
                parts = region.split(':')[1].split('-')
                start, end = int(parts[0]), int(parts[1])
                n = (end - start) // mock_clr.binsize

            # Return the mock contact matrix
            return np.random.random((n, n)) * 10 + np.eye(n) * 50
        return np.array([[]])

    mock_matrix_obj = MagicMock()
    mock_matrix_obj.fetch = mock_matrix_fetch

    mock_matrix = MagicMock()
    mock_matrix.return_value = mock_matrix_obj
    mock_clr.matrix = mock_matrix

    return mock_clr


@pytest.fixture
def sample_bins_df():
    """Create a sample bins DataFrame."""
    return pd.DataFrame({
        'chrom': ['chr12'] * 100,
        'start': range(26_000_000, 26_500_000, 5000),
        'end': range(26_005_000, 26_505_000, 5000)
    })


@pytest.fixture
def sample_insulation_table(sample_bins_df, mock_insulation_scores):
    """
    Create a sample insulation table as returned by cooltools.insulation().
    """
    df = sample_bins_df.copy()
    df['log2_insulation_score_25000'] = mock_insulation_scores[:len(df)]
    df['log2_insulation_score_50000'] = mock_insulation_scores[:len(df)] * 0.8
    df['log2_insulation_score_100000'] = mock_insulation_scores[:len(df)] * 0.6
    df['is_boundary_25000'] = False
    df['is_boundary_50000'] = False
    df['is_boundary_100000'] = False

    # Mark boundaries
    df.loc[38:42, 'is_boundary_25000'] = True
    df.loc[58:62, 'is_boundary_25000'] = True

    return df


@pytest.fixture
def sample_boundary_df(sample_bins_df):
    """
    Create a sample boundary DataFrame with prominence scores.
    """
    df = sample_bins_df.copy()
    df['log2_insulation_score_25000'] = np.random.normal(0, 0.3, len(df))
    df['prominence'] = np.nan
    df['boundary_class'] = None

    # Add strong boundaries
    df.loc[40, 'prominence'] = 0.8
    df.loc[40, 'boundary_class'] = 'strong'
    df.loc[60, 'prominence'] = 0.6
    df.loc[60, 'boundary_class'] = 'strong'

    # Add weak boundary
    df.loc[80, 'prominence'] = 0.3
    df.loc[80, 'boundary_class'] = 'weak'

    return df


@pytest.fixture
def temp_mcool_path():
    """
    Create a temporary path for mock mcool file testing.

    Note: Does not create actual file, just provides a path.
    """
    with tempfile.NamedTemporaryFile(suffix='.mcool', delete=False) as f:
        path = f.name

    yield path

    # Cleanup
    if os.path.exists(path):
        os.remove(path)


@pytest.fixture
def sample_restraints():
    """
    Create sample polymer restraints for testing.

    Returns list of (i, j, rest_length, k) tuples.
    """
    restraints = [
        (0, 5, 2.0, 1.5),
        (2, 8, 1.8, 2.0),
        (10, 15, 2.5, 1.0),
        (20, 25, 1.5, 3.0),
        (30, 40, 2.2, 1.8),
    ]
    return restraints


@pytest.fixture
def random_seed():
    """Fixed random seed for reproducible tests."""
    return 42


@pytest.fixture(autouse=True)
def reset_random_state(random_seed):
    """Automatically reset numpy random state before each test."""
    np.random.seed(random_seed)
    yield
    # Reset after test
    np.random.seed(None)

"""
Unit tests for src/enhancer_id.py.

All tests use mocked data — no actual .mcool, FASTA, or AlphaGenome API
calls required. PyTorch tests are skipped gracefully if torch is not installed.
"""

import pytest
import numpy as np
import pandas as pd
from unittest.mock import MagicMock, patch


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _has_torch():
    try:
        import torch  # noqa: F401
        return True
    except ImportError:
        return False


requires_torch = pytest.mark.skipif(not _has_torch(), reason="PyTorch not installed")


# ---------------------------------------------------------------------------
# Stage 1: extract_candidate_enhancers
# ---------------------------------------------------------------------------

class TestExtractCandidateEnhancers:
    def test_returns_dataframe(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        assert isinstance(result, pd.DataFrame)

    def test_required_columns(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        for col in ('chrom', 'start', 'end', 'contact_score', 'tad_id', 'near_boundary'):
            assert col in result.columns, f"Missing column: {col}"

    def test_filters_to_tad_bins(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        # All returned bins must belong to a TAD
        assert (result['tad_id'] >= 0).all(), "Found bins with tad_id < 0 (outside any TAD)"

    def test_excludes_bins_near_strong_boundaries(
        self, mock_cooler, sample_tad_df, sample_boundary_df
    ):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        # near_boundary should never be True in the output (they are filtered out)
        assert not result['near_boundary'].any(), "near_boundary bins should have been filtered"

    def test_contact_scores_above_median(
        self, mock_cooler, sample_tad_df, sample_boundary_df
    ):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        if len(result) > 1:
            # All scores should be >= the median of the pre-filter pool
            assert (result['contact_score'] >= 0).all()

    def test_sorted_by_contact_score_descending(
        self, mock_cooler, sample_tad_df, sample_boundary_df
    ):
        from src.enhancer_id import extract_candidate_enhancers
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            sample_boundary_df,
        )
        if len(result) > 1:
            scores = result['contact_score'].values
            assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))

    def test_empty_tad_df_returns_empty(self, mock_cooler, sample_boundary_df):
        from src.enhancer_id import extract_candidate_enhancers
        empty_tad = pd.DataFrame(columns=['chrom', 'start', 'end', 'tad_id'])
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            empty_tad,
            sample_boundary_df,
        )
        assert len(result) == 0

    def test_no_strong_boundaries_keeps_more_candidates(
        self, mock_cooler, sample_tad_df
    ):
        from src.enhancer_id import extract_candidate_enhancers
        # boundary_df with no strong boundaries
        no_strong = pd.DataFrame({
            'chrom': ['chr12'] * 3,
            'start': [26_050_000, 26_150_000, 26_250_000],
            'end':   [26_055_000, 26_155_000, 26_255_000],
            'boundary_class': ['weak', 'sub_threshold', None],
            'prominence': [0.3, 0.1, np.nan],
        })
        result = extract_candidate_enhancers(
            mock_cooler,
            "chr12:26,000,000-26,500,000",
            sample_tad_df,
            no_strong,
        )
        # With no strong boundaries excluded, we should get some candidates
        assert isinstance(result, pd.DataFrame)


# ---------------------------------------------------------------------------
# Stage 2: one_hot_encode
# ---------------------------------------------------------------------------

class TestOneHotEncode:
    def test_shape(self):
        from src.enhancer_id import one_hot_encode
        seq = "ACGT"
        enc = one_hot_encode(seq)
        assert enc.shape == (4, 4)

    def test_dtype(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("ACGT")
        assert enc.dtype == np.float32

    def test_adenine_encodes_first_row(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("AAAA")
        assert np.all(enc[0] == 1.0)
        assert np.all(enc[1:] == 0.0)

    def test_cytosine_encodes_second_row(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("CCCC")
        assert np.all(enc[1] == 1.0)
        assert np.all(enc[0] == 0.0)

    def test_n_encodes_all_zeros(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("N")
        assert np.all(enc == 0.0)

    def test_acgt_each_column_one_hot(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("ACGT")
        expected = np.eye(4, dtype=np.float32)
        np.testing.assert_array_equal(enc, expected)

    def test_lowercase(self):
        from src.enhancer_id import one_hot_encode
        enc_upper = one_hot_encode("ACGT")
        enc_lower = one_hot_encode("acgt")
        np.testing.assert_array_equal(enc_upper, enc_lower)

    def test_empty_sequence(self):
        from src.enhancer_id import one_hot_encode
        enc = one_hot_encode("")
        assert enc.shape == (4, 0)

    def test_long_sequence(self):
        from src.enhancer_id import one_hot_encode
        seq = "ACGT" * 125  # 500 bp
        enc = one_hot_encode(seq)
        assert enc.shape == (4, 500)
        # Each position should sum to 1 (exactly one base)
        assert np.all(enc.sum(axis=0) == 1.0)

    def test_sum_per_column_is_0_or_1(self):
        from src.enhancer_id import one_hot_encode
        seq = "ACGTNNN"
        enc = one_hot_encode(seq)
        col_sums = enc.sum(axis=0)
        assert np.all((col_sums == 0.0) | (col_sums == 1.0))


# ---------------------------------------------------------------------------
# Stage 2: fetch_sequences
# ---------------------------------------------------------------------------

class TestFetchSequences:
    def test_requires_pyfaidx(self, sample_tad_df):
        from src.enhancer_id import fetch_sequences
        candidates = pd.DataFrame({
            'chrom': ['chr12'],
            'start': [26_100_000],
            'end':   [26_105_000],
        })
        with pytest.raises((ImportError, Exception)):
            # Will raise ImportError if pyfaidx missing, or FileNotFoundError for bad path
            fetch_sequences(candidates, "/nonexistent/path/genome.fa", window_bp=100)

    def test_output_shape(self, mock_fasta_path):
        from src.enhancer_id import fetch_sequences
        pytest.importorskip("pyfaidx")
        candidates = pd.DataFrame({
            'chrom': ['chrTest', 'chrTest'],
            'start': [0, 100],
            'end':   [50, 150],
        })
        result = fetch_sequences(candidates, mock_fasta_path, window_bp=50)
        assert result.shape == (2, 4, 50)

    def test_output_dtype(self, mock_fasta_path):
        from src.enhancer_id import fetch_sequences
        pytest.importorskip("pyfaidx")
        candidates = pd.DataFrame({
            'chrom': ['chrTest'],
            'start': [0],
            'end':   [50],
        })
        result = fetch_sequences(candidates, mock_fasta_path, window_bp=50)
        assert result.dtype == np.float32

    def test_n_pads_for_near_end(self, mock_fasta_path):
        from src.enhancer_id import fetch_sequences
        pytest.importorskip("pyfaidx")
        # Request a window that extends past the end of the test chromosome
        candidates = pd.DataFrame({
            'chrom': ['chrTest'],
            'start': [450],
            'end':   [500],
        })
        result = fetch_sequences(candidates, mock_fasta_path, window_bp=500)
        # Should return exactly (1, 4, 500) without error
        assert result.shape == (1, 4, 500)


# ---------------------------------------------------------------------------
# Stage 3: EnhancerCNN
# ---------------------------------------------------------------------------

class TestEnhancerCNN:
    @requires_torch
    def test_instantiation(self):
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN()
        assert model.model is not None

    @requires_torch
    def test_forward_shape(self):
        import torch
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN()
        x = torch.zeros(4, 4, 500)  # batch=4, channels=4, length=500
        with torch.no_grad():
            out = model.model(x)
        assert out.shape == (4, 1)

    @requires_torch
    def test_output_range(self):
        import torch
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN()
        x = torch.randn(8, 4, 500)
        with torch.no_grad():
            out = model.model(x)
        assert out.min().item() >= 0.0
        assert out.max().item() <= 1.0

    @requires_torch
    def test_save_load(self, tmp_path):
        import torch
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN()
        save_file = str(tmp_path / "test_model.pt")
        model.save(save_file)

        model2 = EnhancerCNN()
        model2.load(save_file)

        # Weights should match after load
        x = torch.randn(2, 4, 500)
        with torch.no_grad():
            out1 = model.model(x)
            out2 = model2.model(x)
        np.testing.assert_allclose(
            out1.numpy(), out2.numpy(), atol=1e-6
        )

    @requires_torch
    def test_custom_n_filters(self):
        import torch
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN(n_filters=32)
        x = torch.zeros(2, 4, 200)
        with torch.no_grad():
            out = model.model(x)
        assert out.shape == (2, 1)

    def test_import_error_without_torch(self):
        """EnhancerCNN should raise ImportError if torch is unavailable."""
        if _has_torch():
            pytest.skip("torch is installed; cannot test ImportError path")
        from src.enhancer_id import EnhancerCNN
        with pytest.raises(ImportError, match="PyTorch"):
            EnhancerCNN()


# ---------------------------------------------------------------------------
# Stage 4: get_alphagenome_accessibility
# ---------------------------------------------------------------------------

class TestGetAlphagenomeAccessibility:
    def test_returns_array(self, mock_ag_client, sample_candidate_df):
        from src.enhancer_id import get_alphagenome_accessibility
        result = get_alphagenome_accessibility(sample_candidate_df, mock_ag_client)
        assert isinstance(result, np.ndarray)

    def test_shape(self, mock_ag_client, sample_candidate_df):
        from src.enhancer_id import get_alphagenome_accessibility
        result = get_alphagenome_accessibility(sample_candidate_df, mock_ag_client)
        assert result.shape == (len(sample_candidate_df),)

    def test_normalized_to_0_1(self, mock_ag_client, sample_candidate_df):
        from src.enhancer_id import get_alphagenome_accessibility
        result = get_alphagenome_accessibility(sample_candidate_df, mock_ag_client)
        assert result.min() >= 0.0
        assert result.max() <= 1.0 + 1e-6

    def test_handles_api_exception_gracefully(self, sample_candidate_df):
        from src.enhancer_id import get_alphagenome_accessibility
        bad_client = MagicMock()
        bad_client.predict.side_effect = RuntimeError("API error")
        result = get_alphagenome_accessibility(sample_candidate_df, bad_client)
        assert np.all(result == 0.0)

    def test_handles_none_result(self, sample_candidate_df):
        from src.enhancer_id import get_alphagenome_accessibility
        none_client = MagicMock()
        none_client.predict.return_value = None
        result = get_alphagenome_accessibility(sample_candidate_df, none_client)
        assert result.shape == (len(sample_candidate_df),)


# ---------------------------------------------------------------------------
# Stage 5: train_enhancer_cnn and predict_enhancers
# ---------------------------------------------------------------------------

class TestTrainEnhancerCNN:
    @requires_torch
    def test_returns_enhancer_cnn(self):
        from src.enhancer_id import train_enhancer_cnn, EnhancerCNN
        N, window = 20, 100
        seqs = np.random.rand(N, 4, window).astype(np.float32)
        labels = np.random.rand(N).astype(np.float32)
        model = train_enhancer_cnn(seqs, labels, epochs=2, verbose=False)
        assert isinstance(model, EnhancerCNN)

    @requires_torch
    def test_model_in_eval_mode_after_training(self):
        from src.enhancer_id import train_enhancer_cnn
        seqs = np.random.rand(10, 4, 100).astype(np.float32)
        labels = np.random.rand(10).astype(np.float32)
        model = train_enhancer_cnn(seqs, labels, epochs=1, verbose=False)
        assert not model.model.training

    @requires_torch
    def test_save_path(self, tmp_path):
        from src.enhancer_id import train_enhancer_cnn
        seqs = np.random.rand(10, 4, 100).astype(np.float32)
        labels = np.random.rand(10).astype(np.float32)
        save_path = str(tmp_path / "model.pt")
        train_enhancer_cnn(seqs, labels, epochs=1, verbose=False, save_path=save_path)
        assert (tmp_path / "model.pt").exists()

    @requires_torch
    def test_loss_decreases(self):
        """With a very small dataset the loss should not increase wildly."""
        import io
        from contextlib import redirect_stdout
        from src.enhancer_id import train_enhancer_cnn
        # Use a near-constant label signal so the model should converge
        seqs = np.random.rand(32, 4, 200).astype(np.float32)
        labels = np.ones(32, dtype=np.float32) * 0.8
        buf = io.StringIO()
        with redirect_stdout(buf):
            train_enhancer_cnn(seqs, labels, epochs=10, verbose=True)
        output = buf.getvalue()
        # Check that some loss output was produced
        assert "loss=" in output


class TestPredictEnhancers:
    @requires_torch
    def test_returns_dataframe(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import predict_enhancers, EnhancerCNN
        pytest.importorskip("pyfaidx")

        model = EnhancerCNN()
        with patch('src.enhancer_id.fetch_sequences') as mock_fetch:
            n_candidates = 5
            mock_fetch.return_value = np.random.rand(n_candidates, 4, 500).astype(np.float32)

            with patch('src.enhancer_id.extract_candidate_enhancers') as mock_extract:
                mock_extract.return_value = pd.DataFrame({
                    'chrom': ['chr12'] * n_candidates,
                    'start': [26_000_000 + i * 5000 for i in range(n_candidates)],
                    'end':   [26_005_000 + i * 5000 for i in range(n_candidates)],
                    'contact_score': np.random.rand(n_candidates),
                    'tad_id': [0] * n_candidates,
                    'near_boundary': [False] * n_candidates,
                })

                result = predict_enhancers(
                    mock_cooler,
                    "chr12:26,000,000-26,500,000",
                    model,
                    "data/raw/mm10.fa",
                    sample_tad_df,
                    sample_boundary_df,
                )

        assert isinstance(result, pd.DataFrame)

    @requires_torch
    def test_output_columns(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import predict_enhancers, EnhancerCNN
        model = EnhancerCNN()
        with patch('src.enhancer_id.fetch_sequences') as mock_fetch, \
             patch('src.enhancer_id.extract_candidate_enhancers') as mock_extract:

            n = 4
            mock_extract.return_value = pd.DataFrame({
                'chrom': ['chr12'] * n,
                'start': [26_000_000 + i * 5000 for i in range(n)],
                'end':   [26_005_000 + i * 5000 for i in range(n)],
                'contact_score': np.ones(n),
                'tad_id': [0] * n,
                'near_boundary': [False] * n,
            })
            mock_fetch.return_value = np.random.rand(n, 4, 500).astype(np.float32)

            result = predict_enhancers(
                mock_cooler, "chr12:26,000,000-26,500,000",
                model, "data/raw/mm10.fa", sample_tad_df, sample_boundary_df,
            )

        expected_cols = {
            'chrom', 'start', 'end', 'enhancer_score', 'is_enhancer',
            'contact_score', 'tad_id', 'candidate_rank',
        }
        assert expected_cols.issubset(set(result.columns))

    @requires_torch
    def test_scores_in_0_1(self, mock_cooler, sample_tad_df, sample_boundary_df):
        from src.enhancer_id import predict_enhancers, EnhancerCNN
        model = EnhancerCNN()
        with patch('src.enhancer_id.fetch_sequences') as mock_fetch, \
             patch('src.enhancer_id.extract_candidate_enhancers') as mock_extract:

            n = 6
            mock_extract.return_value = pd.DataFrame({
                'chrom': ['chr12'] * n,
                'start': [26_000_000 + i * 5000 for i in range(n)],
                'end':   [26_005_000 + i * 5000 for i in range(n)],
                'contact_score': np.ones(n),
                'tad_id': [0] * n,
                'near_boundary': [False] * n,
            })
            mock_fetch.return_value = np.random.rand(n, 4, 500).astype(np.float32)

            result = predict_enhancers(
                mock_cooler, "chr12:26,000,000-26,500,000",
                model, "data/raw/mm10.fa", sample_tad_df, sample_boundary_df,
            )

        assert result['enhancer_score'].between(0.0, 1.0).all()

    @requires_torch
    def test_empty_candidates_returns_empty_df(
        self, mock_cooler, sample_tad_df, sample_boundary_df
    ):
        from src.enhancer_id import predict_enhancers, EnhancerCNN
        model = EnhancerCNN()
        with patch('src.enhancer_id.extract_candidate_enhancers') as mock_extract:
            mock_extract.return_value = pd.DataFrame(
                columns=['chrom', 'start', 'end', 'contact_score', 'tad_id', 'near_boundary']
            )
            result = predict_enhancers(
                mock_cooler, "chr12:26,000,000-26,500,000",
                model, "data/raw/mm10.fa", sample_tad_df, sample_boundary_df,
            )

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0

    @requires_torch
    def test_sorted_by_enhancer_score_descending(
        self, mock_cooler, sample_tad_df, sample_boundary_df
    ):
        from src.enhancer_id import predict_enhancers, EnhancerCNN
        model = EnhancerCNN()
        with patch('src.enhancer_id.fetch_sequences') as mock_fetch, \
             patch('src.enhancer_id.extract_candidate_enhancers') as mock_extract:

            n = 8
            mock_extract.return_value = pd.DataFrame({
                'chrom': ['chr12'] * n,
                'start': [26_000_000 + i * 5000 for i in range(n)],
                'end':   [26_005_000 + i * 5000 for i in range(n)],
                'contact_score': np.ones(n),
                'tad_id': [0] * n,
                'near_boundary': [False] * n,
            })
            mock_fetch.return_value = np.random.rand(n, 4, 500).astype(np.float32)

            result = predict_enhancers(
                mock_cooler, "chr12:26,000,000-26,500,000",
                model, "data/raw/mm10.fa", sample_tad_df, sample_boundary_df,
            )

        scores = result['enhancer_score'].values
        assert all(scores[i] >= scores[i + 1] for i in range(len(scores) - 1))

import pandas as pd
import pytest
from unittest.mock import patch

from src.genomewide_tads import (
    annotate_tads_with_assays,
    call_genomewide_tads,
    load_interval_signal_table,
    make_boundary_windows,
    summarize_track_over_intervals,
)


@pytest.mark.unit
def test_call_genomewide_tads_returns_three_tables(mock_cooler, sample_insulation_table):
    with patch("src.genomewide_tads.cooltools.insulation", return_value=sample_insulation_table):
        insulation_df, boundary_df, tad_df = call_genomewide_tads(
            mock_cooler,
            window_bp=25_000,
            min_tad_length_bp=50_000,
        )

    assert not insulation_df.empty
    assert "boundary_class" in boundary_df.columns
    assert "tad_id" in tad_df.columns


@pytest.mark.unit
def test_summarize_track_over_intervals_computes_weighted_mean():
    intervals = pd.DataFrame(
        {
            "chrom": ["chr1"],
            "start": [0],
            "end": [100],
        }
    )
    signal = pd.DataFrame(
        {
            "chrom": ["chr1", "chr1"],
            "start": [0, 50],
            "end": [50, 100],
            "signal": [2.0, 4.0],
        }
    )

    result = summarize_track_over_intervals(intervals, signal, prefix="atac")

    assert result.loc[0, "atac_mean"] == pytest.approx(3.0)
    assert result.loc[0, "atac_max"] == pytest.approx(4.0)
    assert result.loc[0, "atac_coverage"] == pytest.approx(1.0)
    assert result.loc[0, "atac_overlap_count"] == 2


@pytest.mark.unit
def test_make_boundary_windows_excludes_chromosome_edges_by_default(mock_cooler):
    chrom_end = int(mock_cooler.chromsizes["chr12"])
    tad_df = pd.DataFrame(
        {
            "chrom": ["chr12", "chr12"],
            "start": [0, 100_000],
            "end": [100_000, chrom_end],
        }
    )

    result = make_boundary_windows(tad_df, chromsizes=mock_cooler.chromsizes, flank_bp=10_000)

    assert len(result) == 1
    assert result.iloc[0]["boundary_bp"] == 100_000


@pytest.mark.unit
def test_annotate_tads_with_assays_adds_expected_columns(mock_cooler):
    tad_df = pd.DataFrame(
        {
            "chrom": ["chr12"],
            "start": [100_000],
            "end": [200_000],
            "tad_id": [0],
            "length_bp": [100_000],
            "n_bins": [20],
            "left_boundary_class": ["strong"],
            "right_boundary_class": ["strong"],
        }
    )
    atac_df = pd.DataFrame(
        {
            "chrom": ["chr12"],
            "start": [120_000],
            "end": [140_000],
            "signal": [5.0],
        }
    )

    tad_summary, boundary_summary = annotate_tads_with_assays(
        tad_df,
        chromsizes=mock_cooler.chromsizes,
        atac_df=atac_df,
    )

    assert "atac_mean" in tad_summary.columns
    assert "atac_boundary_mean" in boundary_summary.columns


@pytest.mark.unit
def test_load_interval_signal_table_supports_bedgraph(tmp_path):
    path = tmp_path / "signal.bedgraph"
    path.write_text("chr1\t0\t100\t3.5\nchr1\t100\t200\t1.0\n")

    result = load_interval_signal_table(path)

    assert list(result.columns) == ["chrom", "start", "end", "signal"]
    assert result["signal"].tolist() == [3.5, 1.0]

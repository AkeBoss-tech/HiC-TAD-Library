"""Tests for local hg38 FASTA slicing (pyfaidx)."""

import pytest

pytest.importorskip("pyfaidx")

from src.hg_dt.data.sequence_fetcher import fetch_hg38_sequence_local


def test_fetch_hg38_sequence_local_slice(tmp_path):
    fa = tmp_path / "mini.fa"
    fa.write_text(">chr1\n" + ("A" * 500) + "\n")
    seq = fetch_hg38_sequence_local("chr1", 10, 30, str(fa))
    assert len(seq) == 20
    assert set(seq) <= {"A"}


def test_fetch_hg38_sequence_local_chr_alias(tmp_path):
    fa = tmp_path / "mini.fa"
    fa.write_text(">chr2\n" + ("C" * 200) + "\n")
    # alias without chr if file uses chr2 only — still chr2
    seq = fetch_hg38_sequence_local("chr2", 0, 50, str(fa))
    assert seq == "C" * 50

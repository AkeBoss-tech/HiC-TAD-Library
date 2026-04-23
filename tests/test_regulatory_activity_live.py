import os

import pytest

from src.hg_dt.models.alphagenome import AlphaGenomeConnector
from src.regulatory_activity import LocusRequest, build_locus_context


pytestmark = [
    pytest.mark.integration,
    pytest.mark.slow,
]


def _has_live_key() -> bool:
    return bool(os.environ.get("ALPHA_GENOME_API_KEY"))


@pytest.mark.skipif(not _has_live_key(), reason="ALPHA_GENOME_API_KEY not set")
def test_live_alphagenome_connector_mouse_rnaseq():
    connector = AlphaGenomeConnector()
    output = connector.predict_interval(
        chrom="chr12",
        start=26_475_712,
        end=27_524_288,
        organism="MOUSE",
        requested_outputs=["RNA_SEQ"],
        cell_type="CL:0000540",
    )
    assert hasattr(output, "rna_seq")
    assert output.rna_seq is not None
    assert len(output.rna_seq.values) > 0


@pytest.mark.skipif(not _has_live_key(), reason="ALPHA_GENOME_API_KEY not set")
def test_live_context_builds_against_real_mouse_assets():
    request = LocusRequest(
        organism="mouse",
        assembly="mm10",
        cell_type="CL:0000540",
        interval="chr12:26000000-28000000",
        fasta_path="data/raw/mm10.fa",
    )
    context = build_locus_context(request)
    assert context.sequence is not None
    assert len(context.sequence) == request.window_bp
    assert context.chrom == "chr12"

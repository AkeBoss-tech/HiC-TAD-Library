import numpy as np
import pandas as pd

from src.regulatory_activity import (
    AlphaGenomePredictionAdapter,
    LocusRequest,
    build_locus_context,
    build_manifest,
    compare_reference_vs_edit,
    load_external_evidence,
    rank_regulatory_activity,
)


class MockTrack:
    def __init__(self, values, resolution=512, metadata=None):
        self.values = np.asarray(values, dtype=float)
        self.resolution = resolution
        self.metadata = metadata if metadata is not None else pd.DataFrame()


class MockOutput:
    def __init__(self, rna_seq=None, cage=None, atac=None, dnase=None, chip_tf=None, contact_maps=None, splice_sites=None, splice_site_usage=None):
        self.rna_seq = rna_seq
        self.cage = cage
        self.atac = atac
        self.dnase = dnase
        self.chip_tf = chip_tf
        self.contact_maps = contact_maps
        self.splice_sites = splice_sites
        self.splice_site_usage = splice_site_usage


class MockConnector:
    def __init__(self, ref_output, edit_output=None):
        self.ref_output = ref_output
        self.edit_output = edit_output if edit_output is not None else ref_output
        self.interval_calls = []
        self.sequence_calls = []

    def predict_interval(self, chrom, start, end, organism="HUMAN", requested_outputs=None, cell_type=None):
        self.interval_calls.append(
            {"chrom": chrom, "start": start, "end": end, "organism": organism, "requested_outputs": requested_outputs, "cell_type": cell_type}
        )
        return self.ref_output

    def predict_sequence(self, sequence, organism="HUMAN", requested_outputs=None, cell_type=None):
        self.sequence_calls.append(
            {"sequence": sequence, "organism": organism, "requested_outputs": requested_outputs, "cell_type": cell_type}
        )
        return self.edit_output if len(self.sequence_calls) > 1 else self.ref_output


def _annotation_df():
    return pd.DataFrame(
        [
            {"feature_type": "gene", "id": "gene1", "Parent": None, "start": 4096, "end": 6656, "strand": "+", "external_name": "GeneA", "chrom": "chr1"},
            {"feature_type": "transcript", "id": "tx1", "Parent": "gene1", "start": 4096, "end": 6656, "strand": "+", "external_name": "GeneA", "chrom": "chr1"},
            {"feature_type": "exon", "id": "ex1", "Parent": "tx1", "start": 4096, "end": 4608, "strand": "+", "external_name": "GeneA", "chrom": "chr1"},
            {"feature_type": "exon", "id": "ex2", "Parent": "tx1", "start": 5632, "end": 6656, "strand": "+", "external_name": "GeneA", "chrom": "chr1"},
        ]
    )


def _make_contact_map(enhancer_boost=1.0):
    n = 16
    mat = np.eye(n) * 3.0
    for i in range(n - 1):
        mat[i, i + 1] = mat[i + 1, i] = 1.2
    mat[4, 8] = mat[8, 4] = 5.0 * enhancer_boost
    mat[3, 8] = mat[8, 3] = 0.8
    mat[10, 8] = mat[8, 10] = 1.0
    return mat


def _mock_output(with_contact=True, with_rna=True, enhancer_boost=1.0):
    n = 16
    atac = np.full(n, 0.05)
    dnase = np.full(n, 0.05)
    cage = np.full(n, 0.03)
    rna = np.full(n, 0.02)
    splice = np.full(n, 0.05)
    ctcf = np.full((n, 2), 0.02)

    # Structural bin
    atac[3] = 0.08
    dnase[3] = 0.04
    ctcf[3, 0] = 1.0
    ctcf[3, 1] = 0.9

    # Enhancer bin
    atac[4] = 1.0 * enhancer_boost
    dnase[4] = 0.8 * enhancer_boost
    cage[4] = 0.15
    rna[4] = 0.12

    # Promoter / gene bin
    cage[8] = 1.0
    rna[8] = 0.95 if with_rna else 0.0
    splice[8] = 0.7

    metadata = pd.DataFrame({"transcription_factor": ["CTCF", "OTHER"]})
    return MockOutput(
        rna_seq=MockTrack(rna if with_rna else np.zeros(n), resolution=512),
        cage=MockTrack(cage, resolution=512),
        atac=MockTrack(atac, resolution=512),
        dnase=MockTrack(dnase, resolution=512),
        chip_tf=MockTrack(ctcf, resolution=512, metadata=metadata),
        contact_maps=MockTrack(_make_contact_map(enhancer_boost=enhancer_boost), resolution=512) if with_contact else None,
        splice_sites=MockTrack(splice, resolution=512),
        splice_site_usage=MockTrack(splice, resolution=512),
    )


def _build_context(sequence=None, organism="mouse", cell_type="neuron"):
    request = LocusRequest(
        organism=organism,
        assembly="mm10" if organism == "mouse" else "hg38",
        cell_type=cell_type,
        interval="chr1:2048-6144",
        sequence=sequence if sequence is not None else ("A" * 8192),
        anchor_interval="chr1:0-8192",
        window_bp=8192,
        ranking_resolution_bp=512,
    )
    return build_locus_context(request, annotation_df=_annotation_df())


def test_builds_locus_context_from_interval_with_transcripts():
    context = _build_context()
    assert context.chrom == "chr1"
    assert len(context.bins) == 16
    assert set(context.transcripts["transcript_id"]) == {"tx1"}
    assert set(context.genes["gene_id"]) == {"gene1"}


def test_ranks_promoter_and_enhancer_above_inert_bins():
    context = _build_context()
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output()))
    prediction = adapter.predict_reference(context)
    bins, links, genes, summary = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    promoter_row = bins[bins["bin_index"] == 8].iloc[0]
    enhancer_row = bins[bins["bin_index"] == 4].iloc[0]
    inert_row = bins[bins["bin_index"] == 1].iloc[0]

    assert promoter_row["element_type"] == "promoter"
    assert enhancer_row["activity_score"] > inert_row["activity_score"]
    assert summary["n_links"] >= 1
    assert genes[0].gene_symbol == "GeneA"


def test_structural_bins_rank_below_true_enhancers_without_sequence_support():
    context = _build_context()
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output()))
    prediction = adapter.predict_reference(context)
    bins, _, _, _ = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    structural = bins[bins["bin_index"] == 3].iloc[0]
    enhancer = bins[bins["bin_index"] == 4].iloc[0]
    assert structural["element_type"] == "structural"
    assert structural["activity_score"] < enhancer["activity_score"]


def test_links_enhancer_to_correct_promoter_using_contact_support():
    context = _build_context()
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output()))
    prediction = adapter.predict_reference(context)
    bins, links, _, _ = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    best_link = links.iloc[0]
    enhancer_row = bins[bins["bin_index"] == 4].iloc[0]
    assert int(best_link["enhancer_bin"]) == 4
    assert int(best_link["promoter_bin"]) == 8
    assert enhancer_row["linked_gene_symbols"] == ["GeneA"]


def test_stable_when_rna_present_but_contact_missing():
    context = _build_context()
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output(with_contact=False)))
    prediction = adapter.predict_reference(context)
    bins, links, _, _ = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    assert "contact_maps" in prediction.provenance["missing_modalities"]
    assert not bins["activity_score"].isna().any()
    assert links is not None


def test_stable_when_contact_present_but_external_assay_missing():
    context = _build_context()
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output()))
    prediction = adapter.predict_reference(context)
    bins, _, _, _ = rank_regulatory_activity(context, prediction, {})
    assert not bins["cell_type_support_score"].isna().any()
    assert np.all(bins["cell_type_support_score"].values >= 0.0)


def test_perturbation_reports_loss_after_synthetic_enhancer_deletion():
    context = _build_context()
    ref_output = _mock_output()
    edit_output = _mock_output(enhancer_boost=0.2)
    adapter = AlphaGenomePredictionAdapter(MockConnector(ref_output, edit_output))
    prediction = adapter.predict_reference(context)
    ref_bins, _, _, _ = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    edited_bins, _, edited_genes, _, diff_summary, _ = compare_reference_vs_edit(
        context,
        adapter,
        ref_bins,
        {"mode": "deletion", "start": 2048, "end": 2560},
        {"observed_contact": np.zeros(len(context.bins))},
    )

    enhancer_delta = edited_bins[edited_bins["bin_index"] == 4].iloc[0]["delta_activity_score"]
    assert enhancer_delta < 0.0
    assert 4 in diff_summary["lost_bins"]
    assert edited_genes[0].linked_element_score <= 1.0


def test_perturbation_reports_gain_after_synthetic_enhancer_insertion():
    context = _build_context()
    ref_output = _mock_output(enhancer_boost=0.2)
    edit_output = _mock_output(enhancer_boost=1.3)
    adapter = AlphaGenomePredictionAdapter(MockConnector(ref_output, edit_output))
    prediction = adapter.predict_reference(context)
    ref_bins, _, _, _ = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})

    edited_bins, _, _, _, diff_summary, _ = compare_reference_vs_edit(
        context,
        adapter,
        ref_bins,
        {"mode": "insertion", "start": 2048, "end": 2048, "sequence": "A" * 64},
        {"observed_contact": np.zeros(len(context.bins))},
    )

    enhancer_delta = edited_bins[edited_bins["bin_index"] == 4].iloc[0]["delta_activity_score"]
    assert enhancer_delta > 0.0
    assert 4 in diff_summary["gained_bins"]


def test_mouse_and_human_requests_route_to_correct_organisms():
    mouse_context = _build_context(organism="mouse")
    human_context = _build_context(organism="human")
    mouse_connector = MockConnector(_mock_output())
    human_connector = MockConnector(_mock_output())

    AlphaGenomePredictionAdapter(mouse_connector).predict_reference(mouse_context)
    AlphaGenomePredictionAdapter(human_connector).predict_reference(human_context)

    assert mouse_connector.sequence_calls[0]["organism"] == "MOUSE"
    assert human_connector.sequence_calls[0]["organism"] == "HUMAN"


def test_manifest_records_missing_modalities_and_biosample_fallbacks():
    context = _build_context(cell_type="totally_unknown_celltype")
    adapter = AlphaGenomePredictionAdapter(MockConnector(_mock_output(with_contact=False)))
    prediction = adapter.predict_reference(context)
    bins, links, genes, summary = rank_regulatory_activity(context, prediction, {"observed_contact": np.zeros(len(context.bins))})
    _, external_provenance = load_external_evidence(context)
    manifest = build_manifest(context, bins, links, genes, prediction, external_provenance, summary)

    assert "contact_maps" in manifest["model_provenance"]["predicted"]["missing_modalities"]
    assert "external" in manifest["model_provenance"]
    assert isinstance(manifest["model_provenance"]["external"]["biosample_fallbacks"], list)

import os
import pytest
import time
import pandas as pd
from src.hg_dt.data.builder import ReferenceContextBuilder

@pytest.fixture
def mock_data_dir(tmpdir):
    data_dir = tmpdir.mkdir("mock_data")

    # Create mock FASTA
    # 2 Mb sequence of 'A's, 'C's, 'G's, 'T's for chr17
    # We will simulate BRCA1 locus at chr17:43000000 (0-based)
    # Actually, let's make a smaller mock fasta that covers a 2Mb range
    # e.g., 0 to 2,000,000 for simplicity and map our "chr17" coordinates to 0-2Mb.

    fasta_path = os.path.join(data_dir, "test.fa")
    seq_length = 2000000
    with open(fasta_path, "w") as f:
        f.write(">chr17\n")
        # Write chunks of 80 bases
        full_seq = "ACGT" * (seq_length // 4)
        for i in range(0, seq_length, 80):
            f.write(full_seq[i:i+80] + "\n")

    # Create mock GENCODE GTF
    gtf_path = os.path.join(data_dir, "test.gtf")
    with open(gtf_path, "w") as f:
        # Chromosome, Source, Feature, Start, End, Score, Strand, Frame, Attribute
        f.write("chr17\tHAVANA\tgene\t1000000\t1010000\t.\t+\t.\tgene_id \"ENSG00000012048.24\"; gene_name \"BRCA1\";\n")
        f.write("chr17\tHAVANA\texon\t1000000\t1000500\t.\t+\t.\tgene_id \"ENSG00000012048.24\"; gene_name \"BRCA1\";\n")

    # Create mock ENCODE BED
    bed_path = os.path.join(data_dir, "test.bed")
    with open(bed_path, "w") as f:
        # chr, start, end, name
        f.write("chr17\t995000\t1000000\tEH38E2776519\t255\t.\t995000\t1000000\t255,0,0\n")

    return fasta_path, gtf_path, bed_path

def test_reference_context_builder(mock_data_dir):
    fasta_path, gtf_path, bed_path = mock_data_dir
    builder = ReferenceContextBuilder(fasta_path, gtf_path, bed_path)

    # Test execution speed requirement (< 200ms)
    start_time = time.time()

    # Query a 1Mb window around chr17:1000000
    # Let's say there is a 5 kb enhancer deletion from 995000 to 1000000
    edit_start = 995000
    edit_end = 1000000

    context = builder.get_context(
        chrom="chr17",
        start=edit_start,
        end=edit_end,
        edit_type="deletion",
        window_size=1000000
    )

    exec_time = time.time() - start_time
    assert exec_time < 0.2, f"Execution took too long: {exec_time}s"

    # Check reference and mutant sequences
    ref_seq = context['ref_seq']
    mut_seq = context['mut_seq']

    assert len(ref_seq) == 1000000
    assert len(mut_seq) == 1000000 - 5000

    # Verify annotations
    annotations = context['annotations']
    genes_df = annotations['genes']
    ccres_df = annotations['cCREs']

    assert not genes_df.empty, "Should identify exons/genes"

    # GTF reader in pyranges parses attributes into separate columns
    if 'gene_name' in genes_df.columns:
        assert "BRCA1" in genes_df['gene_name'].values
    else:
        assert "BRCA1" in str(genes_df.to_dict())

    assert not ccres_df.empty, "Should identify enhancers"
    # Ensure our mock enhancer is in there
    # BED file columns might be automatically assigned if headerless
    # but the start coordinate should be there
    assert 995000 in ccres_df['Start'].values

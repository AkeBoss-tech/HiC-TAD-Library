#!/usr/bin/env python3
"""
Check what output types AlphaGenome has available for mouse cell types.
"""
import os
from dotenv import load_dotenv
from alphagenome.data import genome
from alphagenome.models import dna_client

load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))

# Test region
interval = genome.Interval(chromosome='chr13', start=83_218_179, end=84_266_755)

# Cell types
cell_types = {
    'CL:0000207': 'Olfactory receptor cell',
    'EFO:0004038': 'Mouse embryonic stem cell',
}

print("Checking available data for mouse cell types...\n")

for cell_type, name in cell_types.items():
    print(f"\n{'='*70}")
    print(f"{cell_type} — {name}")
    print(f"{'='*70}")

    # Request all output types
    result = dna_model.predict_interval(
        interval=interval,
        requested_outputs={
            dna_client.OutputType.CONTACT_MAPS,
            dna_client.OutputType.CHIP_TF,
            dna_client.OutputType.ATAC,
            dna_client.OutputType.RNA_SEQ,
            dna_client.OutputType.DNASE,
            dna_client.OutputType.CHIP_HISTONE,
        },
        ontology_terms=[cell_type],
        organism=dna_client.Organism.MUS_MUSCULUS,
    )

    # Check what we got
    print("\nData availability:")

    if result.contact_maps is not None:
        print(f"  ✓ Contact Maps: {result.contact_maps.values.shape}")
        print(f"    Cell types: {len(result.contact_maps.ontology_terms)}")
    else:
        print(f"  ✗ Contact Maps: Not available")

    if result.chip_tf is not None:
        print(f"  ✓ TF Binding (CHIP_TF): {result.chip_tf.values.shape}")
        if result.chip_tf.ontology_terms is not None:
            print(f"    Cell types: {len(result.chip_tf.ontology_terms)}")
        if result.chip_tf.values.shape[1] > 0:
            print(f"    Non-zero values: {(result.chip_tf.values > 0).sum()}")
        else:
            print(f"    ⚠ No TF tracks available (shape shows 0 tracks)")
    else:
        print(f"  ✗ TF Binding (CHIP_TF): Not available")

    if result.atac is not None:
        print(f"  ✓ ATAC-seq: {result.atac.values.shape}")
        if result.atac.ontology_terms is not None:
            print(f"    Cell types: {len(result.atac.ontology_terms)}")
        if result.atac.values.shape[1] > 0:
            print(f"    Non-zero values: {(result.atac.values > 0).sum()}")
        else:
            print(f"    ⚠ No ATAC tracks available (shape shows 0 tracks)")
    else:
        print(f"  ✗ ATAC-seq: Not available")

    if result.rna_seq is not None:
        print(f"  ✓ RNA-seq: {result.rna_seq.values.shape}")
        if result.rna_seq.ontology_terms is not None:
            print(f"    Cell types: {len(result.rna_seq.ontology_terms)}")
        if result.rna_seq.values.shape[1] == 0:
            print("    ⚠ No RNA-seq tracks available")
    else:
        print("  ✗ RNA-seq: Not available")

    if result.dnase is not None:
        print(f"  ✓ DNase: {result.dnase.values.shape}")
        if result.dnase.ontology_terms is not None:
            print(f"    Cell types: {len(result.dnase.ontology_terms)}")
        if result.dnase.values.shape[1] == 0:
            print("    ⚠ No DNase tracks available")
    else:
        print("  ✗ DNase: Not available")

    if result.chip_histone is not None:
        print(f"  ✓ Histone ChIP: {result.chip_histone.values.shape}")
        if result.chip_histone.ontology_terms is not None:
            print(f"    Cell types: {len(result.chip_histone.ontology_terms)}")
        if result.chip_histone.values.shape[1] == 0:
            print("    ⚠ No histone ChIP tracks available")
    else:
        print("  ✗ Histone ChIP: Not available")

print("\n" + "="*70)
print("Summary:")
print("="*70)
print("AlphaGenome mouse data is limited compared to human.")
print("Most mouse tracks are only available for a few cell types.")
print("Consider using human cell types (e.g., K562, H1-hESC) for full multi-modal analysis.")

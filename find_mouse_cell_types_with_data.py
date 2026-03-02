#!/usr/bin/env python3
"""
Scan through common mouse cell types to find which have TF binding and ATAC data.
"""
import os
from dotenv import load_dotenv
from alphagenome.data import genome
from alphagenome.models import dna_client

load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))

# Test region
interval = genome.Interval(chromosome='chr13', start=83_218_179, end=84_266_755)

# Common mouse cell type ontology terms
# Source: Cell Ontology and EFO
mouse_cell_types = {
    # Stem cells
    'EFO:0004038': 'Mouse embryonic stem cell (mESC)',
    'CL:0002322': 'Embryonic stem cell',

    # Neural
    'CL:0000540': 'Neuron',
    'CL:0000125': 'Glial cell',
    'CL:0000617': 'GABAergic neuron',
    'CL:0000700': 'Dopaminergic neuron',

    # Immune
    'CL:0000542': 'Lymphocyte',
    'CL:0000623': 'Natural killer cell',
    'CL:0000624': 'CD4-positive T cell',
    'CL:0000625': 'CD8-positive T cell',
    'CL:0000235': 'Macrophage',
    'CL:0000451': 'Dendritic cell',

    # Other
    'CL:0000182': 'Hepatocyte',
    'CL:0000746': 'Cardiac muscle cell',
    'CL:0002494': 'Cardiomyocyte',
    'CL:0000187': 'Muscle cell',
    'CL:0000057': 'Fibroblast',
    'CL:0000115': 'Endothelial cell',
    'CL:0000670': 'Primordial germ cell',

    # Sensory
    'CL:0000207': 'Olfactory receptor cell',
}

print("Scanning mouse cell types for available data...\n")
print(f"{'Cell Type':<50} {'Contact':<8} {'TF':<8} {'ATAC':<8}")
print("="*74)

available_multimodal = []

for cell_type, name in sorted(mouse_cell_types.items()):
    try:
        result = dna_model.predict_interval(
            interval=interval,
            requested_outputs={
                dna_client.OutputType.CONTACT_MAPS,
                dna_client.OutputType.CHIP_TF,
                dna_client.OutputType.ATAC,
            },
            ontology_terms=[cell_type],
            organism=dna_client.Organism.MUS_MUSCULUS,
        )

        has_contact = result.contact_maps.values.shape[-1] > 0
        has_tf = result.chip_tf.values.shape[1] > 0
        has_atac = result.atac.values.shape[1] > 0

        contact_str = "✓" if has_contact else "✗"
        tf_str = "✓" if has_tf else "✗"
        atac_str = "✓" if has_atac else "✗"

        print(f"{name:<50} {contact_str:<8} {tf_str:<8} {atac_str:<8}")

        if has_tf and has_atac:
            available_multimodal.append((cell_type, name))

    except Exception as e:
        print(f"{name:<50} ERROR: {str(e)[:20]}")

print("\n" + "="*74)
print("\nCell types with BOTH TF binding AND ATAC data:")
print("="*74)

if available_multimodal:
    for cell_type, name in available_multimodal:
        print(f"  {cell_type}: {name}")
else:
    print("  ⚠ No mouse cell types found with both TF and ATAC data")
    print("  Consider using human cell types (HOMO_SAPIENS) instead")

print("\n" + "="*74)
print("Next steps:")
print("="*74)
print("1. If cell types found: Update CELL_TYPES in run_mouse_deletion.py")
print("2. If no cell types found: Switch to human organism with K562, H1-hESC, etc.")

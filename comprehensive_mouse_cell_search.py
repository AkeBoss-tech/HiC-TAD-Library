#!/usr/bin/env python3
"""
Comprehensive search for mouse cell types with multi-modal data.
Tests a wide range of ontology terms to find ANY cell type with TF/ATAC data.
"""
import os
from dotenv import load_dotenv
from alphagenome.data import genome
from alphagenome.models import dna_client

load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))

# Test region
interval = genome.Interval(chromosome='chr13', start=83_218_179, end=84_266_755)

# Comprehensive list of mouse cell type ontology terms
# Includes EFO, CL, UBERON, and other ontologies
mouse_cell_types = {
    # ── Stem & Progenitor Cells ────────────────────────────────────────────────
    'EFO:0004038': 'Mouse embryonic stem cell (mESC)',
    'CL:0002322': 'Embryonic stem cell',
    'CL:0000034': 'Stem cell',
    'CL:0002248': 'Pluripotent stem cell',
    'CL:0000048': 'Multi-fate stem cell',

    # ── Neural & Brain ──────────────────────────────────────────────────────────
    'CL:0000540': 'Neuron',
    'CL:0000125': 'Glial cell',
    'CL:0000617': 'GABAergic neuron',
    'CL:0000700': 'Dopaminergic neuron',
    'CL:0000746': 'Cardiac muscle cell',
    'UBERON:0000955': 'Brain',
    'UBERON:0001851': 'Cortex',
    'UBERON:0002421': 'Hippocampus',
    'UBERON:0002037': 'Cerebellum',

    # ── Immune System ───────────────────────────────────────────────────────────
    'CL:0000542': 'Lymphocyte',
    'CL:0000623': 'Natural killer cell',
    'CL:0000624': 'CD4-positive T cell',
    'CL:0000625': 'CD8-positive T cell',
    'CL:0000235': 'Macrophage',
    'CL:0000451': 'Dendritic cell',
    'CL:0000236': 'B cell',
    'CL:0000815': 'Regulatory T cell',
    'CL:0000066': 'Epithelial cell',
    'UBERON:0002371': 'Bone marrow',
    'UBERON:0002106': 'Spleen',
    'UBERON:0002370': 'Thymus',

    # ── Liver & Metabolic ───────────────────────────────────────────────────────
    'CL:0000182': 'Hepatocyte',
    'UBERON:0002107': 'Liver',

    # ── Heart & Muscle ──────────────────────────────────────────────────────────
    'CL:0002494': 'Cardiomyocyte',
    'CL:0000187': 'Muscle cell',
    'CL:0000359': 'Vascular associated smooth muscle cell',
    'UBERON:0000948': 'Heart',
    'UBERON:0001134': 'Skeletal muscle tissue',

    # ── Connective Tissue ───────────────────────────────────────────────────────
    'CL:0000057': 'Fibroblast',
    'CL:0000115': 'Endothelial cell',
    'CL:0000136': 'Fat cell',

    # ── Reproductive & Development ──────────────────────────────────────────────
    'CL:0000670': 'Primordial germ cell',
    'UBERON:0000922': 'Embryo',
    'UBERON:0007023': 'Adult organism',

    # ── Sensory & Specialized ───────────────────────────────────────────────────
    'CL:0000207': 'Olfactory receptor cell',
    'CL:0000210': 'Photoreceptor cell',

    # ── Kidney & Excretory ──────────────────────────────────────────────────────
    'UBERON:0002113': 'Kidney',
    'CL:0000653': 'Podocyte',

    # ── Lung & Respiratory ──────────────────────────────────────────────────────
    'UBERON:0002048': 'Lung',
    'CL:0002062': 'Type II pneumocyte',

    # ── Intestine & Digestive ───────────────────────────────────────────────────
    'UBERON:0000160': 'Intestine',
    'CL:0000584': 'Enterocyte',

    # ── Development Stages ──────────────────────────────────────────────────────
    'UBERON:0000068': 'Embryonic structure',
    'UBERON:0007220': 'E10.5 embryo',
    'UBERON:0007221': 'E11.5 embryo',
    'UBERON:0007222': 'E12.5 embryo',
    'UBERON:0007223': 'E13.5 embryo',
    'UBERON:0007224': 'E14.5 embryo',
    'UBERON:0007225': 'E15.5 embryo',
}

print("="*80)
print("COMPREHENSIVE MOUSE CELL TYPE SCAN FOR MULTI-MODAL DATA")
print("="*80)
print(f"Testing {len(mouse_cell_types)} cell types/tissues/developmental stages\n")

# Track results
contact_only = []
has_tf = []
has_atac = []
has_both = []
has_either = []
has_rna = []
errors = []

print(f"{'ID':<20} {'Name':<40} {'Contact':<8} {'TF':<6} {'ATAC':<6} {'RNA':<6}")
print("="*80)

for cell_id, name in sorted(mouse_cell_types.items(), key=lambda x: x[1]):
    try:
        result = dna_model.predict_interval(
            interval=interval,
            requested_outputs={
                dna_client.OutputType.CONTACT_MAPS,
                dna_client.OutputType.CHIP_TF,
                dna_client.OutputType.ATAC,
                dna_client.OutputType.RNA_SEQ,
            },
            ontology_terms=[cell_id],
            organism=dna_client.Organism.MUS_MUSCULUS,
        )

        has_contact = result.contact_maps.values.shape[-1] > 0
        has_tf_data = result.chip_tf.values.shape[1] > 0
        has_atac_data = result.atac.values.shape[1] > 0
        has_rna_data = result.rna_seq.values.shape[1] > 0

        c_str = "✓" if has_contact else "✗"
        t_str = "✓" if has_tf_data else "✗"
        a_str = "✓" if has_atac_data else "✗"
        r_str = "✓" if has_rna_data else "✗"

        print(f"{cell_id:<20} {name[:38]:<40} {c_str:<8} {t_str:<6} {a_str:<6} {r_str:<6}")

        # Categorize
        if has_tf_data and has_atac_data:
            has_both.append((cell_id, name))
        elif has_tf_data or has_atac_data:
            has_either.append((cell_id, name, 'TF' if has_tf_data else 'ATAC'))
        elif has_contact:
            contact_only.append((cell_id, name))

        if has_tf_data:
            has_tf.append((cell_id, name))
        if has_atac_data:
            has_atac.append((cell_id, name))
        if has_rna_data:
            has_rna.append((cell_id, name))

    except Exception as e:
        errors.append((cell_id, name, str(e)[:30]))
        print(f"{cell_id:<20} {name[:38]:<40} ERROR: {str(e)[:20]}")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

print(f"\n✓ Cell types with BOTH TF and ATAC: {len(has_both)}")
if has_both:
    for cell_id, name in has_both:
        print(f"  • {cell_id}: {name}")
else:
    print("  (None found)")

print(f"\n✓ Cell types with EITHER TF or ATAC: {len(has_either)}")
if has_either:
    for cell_id, name, modality in has_either:
        print(f"  • {cell_id}: {name} ({modality} available)")
else:
    print("  (None found)")

print(f"\n✓ Cell types with TF binding: {len(has_tf)}")
if has_tf:
    for cell_id, name in has_tf[:10]:  # Show first 10
        print(f"  • {cell_id}: {name}")
    if len(has_tf) > 10:
        print(f"  ... and {len(has_tf) - 10} more")

print(f"\n✓ Cell types with ATAC-seq: {len(has_atac)}")
if has_atac:
    for cell_id, name in has_atac[:10]:
        print(f"  • {cell_id}: {name}")
    if len(has_atac) > 10:
        print(f"  ... and {len(has_atac) - 10} more")

print(f"\n✓ Cell types with RNA-seq: {len(has_rna)}")
if has_rna:
    for cell_id, name in has_rna[:10]:
        print(f"  • {cell_id}: {name}")
    if len(has_rna) > 10:
        print(f"  ... and {len(has_rna) - 10} more")

print(f"\n✓ Cell types with contact maps only: {len(contact_only)}")

print(f"\n✗ Errors encountered: {len(errors)}")
if errors and len(errors) <= 5:
    for cell_id, name, err in errors:
        print(f"  • {cell_id}: {err}")

print("\n" + "="*80)
print("RECOMMENDATION")
print("="*80)

if has_both:
    print("✓ Use these cell types for FULL multi-modal analysis (Contact + TF + ATAC):")
    for cell_id, name in has_both:
        print(f"\n  CELL_TYPES = {{")
        print(f"      '{cell_id}': '{name}',")
        print(f"  }}")
elif has_either:
    print("⚠ No cell types with BOTH TF and ATAC, but some have one or the other:")
    print("\nCell types with TF binding:")
    for cell_id, name in has_tf[:3]:
        print(f"  • '{cell_id}': '{name}'")
    print("\nCell types with ATAC-seq:")
    for cell_id, name in has_atac[:3]:
        print(f"  • '{cell_id}': '{name}'")
    print("\nYou could run separate analyses for TF and ATAC.")
else:
    print("✗ No mouse cell types found with TF binding or ATAC-seq data.")
    print("  Recommendation: Use human cell types (HOMO_SAPIENS) for multi-modal analysis.")
    print("\n  Example human cell types with full data:")
    print("    • EFO:0002067: K562 (erythroleukemia)")
    print("    • EFO:0003042: H1-hESC (embryonic stem cells)")
    print("    • EFO:0001200: HepG2 (liver)")
    print("    • EFO:0002784: GM12878 (lymphoblastoid)")

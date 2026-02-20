#!/usr/bin/env python3
"""
Find mouse cell types/tissues with CONTACT MAPS + TF + ATAC (all three).
"""
import os
from dotenv import load_dotenv
from alphagenome.data import genome
from alphagenome.models import dna_client

load_dotenv()
dna_model = dna_client.create(os.getenv('ALPHA_GENOME_API_KEY'))

interval = genome.Interval(chromosome='chr13', start=83_218_179, end=84_266_755)

# Results from previous scan
candidates = {
    # Tissues with TF + ATAC (need to check contact maps)
    'UBERON:0000948': 'Heart',
    'UBERON:0000160': 'Intestine',
    'UBERON:0002113': 'Kidney',
    'UBERON:0002107': 'Liver',
    'UBERON:0002048': 'Lung',

    # Tissues with contact maps (need to check TF/ATAC)
    'EFO:0004038': 'Mouse embryonic stem cell',
    'CL:0000207': 'Olfactory receptor cell',

    # Additional brain/tissue terms
    'UBERON:0000955': 'Brain',
    'UBERON:0002371': 'Bone marrow',
    'UBERON:0002037': 'Cerebellum',
    'UBERON:0002106': 'Spleen',
}

print("="*80)
print("SEARCHING FOR MOUSE TISSUES WITH CONTACT MAPS + TF + ATAC")
print("="*80)
print("\nTesting candidates...\n")

found_triple = []

for cell_id, name in sorted(candidates.items()):
    try:
        result = dna_model.predict_interval(
            interval=interval,
            requested_outputs={
                dna_client.OutputType.CONTACT_MAPS,
                dna_client.OutputType.CHIP_TF,
                dna_client.OutputType.ATAC,
            },
            ontology_terms=[cell_id],
            organism=dna_client.Organism.MUS_MUSCULUS,
        )

        has_contact = result.contact_maps.values.shape[-1] > 0
        has_tf = result.chip_tf.values.shape[1] > 0
        has_atac = result.atac.values.shape[1] > 0

        contact_count = result.contact_maps.values.shape[-1] if has_contact else 0
        tf_count = result.chip_tf.values.shape[1] if has_tf else 0
        atac_count = result.atac.values.shape[1] if has_atac else 0

        status = "✓ ALL THREE" if (has_contact and has_tf and has_atac) else "✗ Missing"

        print(f"{cell_id:<20} {name:<30} {status}")
        print(f"  Contact: {contact_count:>3} tracks | TF: {tf_count:>3} tracks | ATAC: {atac_count:>3} tracks")

        if has_contact and has_tf and has_atac:
            found_triple.append((cell_id, name, contact_count, tf_count, atac_count))

    except Exception as e:
        print(f"{cell_id:<20} {name:<30} ERROR: {str(e)[:40]}")

print("\n" + "="*80)
print("RESULTS")
print("="*80)

if found_triple:
    print(f"\n✓ Found {len(found_triple)} mouse tissue(s) with ALL THREE modalities:\n")
    for cell_id, name, c_cnt, t_cnt, a_cnt in found_triple:
        print(f"  {cell_id}: {name}")
        print(f"    - Contact maps: {c_cnt} tracks")
        print(f"    - TF binding: {t_cnt} tracks")
        print(f"    - ATAC-seq: {a_cnt} tracks")
        print()

    print("RECOMMENDED CONFIGURATION:")
    print("="*80)
    for cell_id, name, _, _, _ in found_triple:
        print(f"\nCELL_TYPES = {{")
        print(f"    '{cell_id}': '{name}',")
        print(f"}}")
        print(f"\n# Use this in run_mouse_deletion.py for full multi-modal analysis!")

else:
    print("\n✗ No mouse tissues found with ALL THREE: contact maps + TF + ATAC")
    print("\nCurrent situation:")
    print("  • Mouse ESC & Olfactory: Have contact maps, NO TF/ATAC")
    print("  • Heart/Liver/Lung/etc: Have TF + ATAC, NO contact maps")
    print("\nOPTIONS:")
    print("  1. Use Mouse ESC/Olfactory for contact map analysis (current approach)")
    print("  2. Switch to human (HOMO_SAPIENS) for full multi-modal with contact maps")
    print("\nHuman cell types with ALL THREE:")
    print("  • EFO:0002067: K562 (erythroleukemia)")
    print("  • EFO:0003042: H1-hESC (human embryonic stem cells)")

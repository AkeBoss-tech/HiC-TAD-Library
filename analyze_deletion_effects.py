#!/usr/bin/env python3
"""
Analyze DNA deletion effects on chromatin contact maps using AlphaGenome.

This script demonstrates how structural variants (deletions) affect 3D genome
organization by:
1. Predicting contact maps for wild-type mouse genomic regions
2. Simulating deletions (e.g., CTCF sites, TAD boundaries, regulatory elements)
3. Predicting contact maps for deleted sequences
4. Visualizing before/after comparison

The deletions can reveal:
- TAD boundary disruption effects
- CTCF-dependent loop anchoring
- Regulatory element contributions to chromatin architecture
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from dotenv import load_dotenv

from alphagenome.data import genome
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# Load environment variables
load_dotenv()
API_KEY = os.getenv('ALPHA_GENOME_API_KEY')

# Initialize AlphaGenome model
print("Initializing AlphaGenome model...")
dna_model = dna_client.create(API_KEY)


def visualize_deletion_effect(
    interval_str,
    deletion_start,
    deletion_end,
    cell_type=None,  # Use all available cell types for mouse
    organism=dna_client.Organism.MUS_MUSCULUS,
    title_prefix="TAD Boundary Deletion"
):
    """
    Visualize the effect of a DNA deletion on contact maps.

    Args:
        interval_str: Genomic interval (e.g., "chr12:26000000-27000000")
        deletion_start: Start position of deletion within interval
        deletion_end: End position of deletion within interval
        cell_type: Cell type ontology term (default: K562)
        organism: Organism (default: mouse)
        title_prefix: Prefix for plot titles
    """
    print(f"\n{'='*70}")
    print(f"Analyzing deletion effect in {interval_str}")
    print(f"Deletion region: {deletion_start}-{deletion_end}")
    print(f"Deletion size: {deletion_end - deletion_start:,} bp")
    print(f"{'='*70}\n")

    # Parse interval
    interval = genome.Interval.from_str(interval_str)

    # Resize to model's supported length (2^20 = 1,048,576 bp)
    resized_interval = interval.resize(2**20)
    print(f"Resized interval for model: {resized_interval}")

    # Create deletion variant
    # AlphaGenome uses VCF-style variant notation where:
    # - For deletions, ref is the deleted sequence + flanking base
    # - Alt is just the flanking base
    # For simplicity, we'll create a large deletion variant
    deletion_variant = genome.Variant(
        chromosome=interval.chromosome,
        position=deletion_start - 1,  # 0-based, one before deletion
        reference_bases='N' * (deletion_end - deletion_start + 1),  # Placeholder
        alternate_bases='N',  # Single base (represents deletion)
    )

    print(f"\n1. Predicting WILD-TYPE contact map...")
    try:
        # Predict wild-type contact map
        # Try without ontology filter first to see what's available
        wt_output = dna_model.predict_interval(
            interval=resized_interval,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            organism=organism,
            ontology_terms=None  # Get all available cell types
        )
        print("   ✓ Wild-type prediction complete")

        # Extract contact matrix
        wt_contact_map = wt_output.contact_maps
        print(f"   Contact map shape: {wt_contact_map.values.shape}")
        print(f"   Available cell types: {len(wt_contact_map.ontology_terms)}")
        if len(wt_contact_map.ontology_terms) > 0:
            print(f"   Cell type IDs: {[term.ontology_curie for term in wt_contact_map.ontology_terms[:5]]}")

        # Check if we got any data
        if wt_contact_map.values.shape[-1] == 0:
            print(f"   ⚠ No contact map data available for specified cell type")
            print(f"   Trying without cell type filter to find available data...")
            # Retry without cell type filter
            wt_output = dna_model.predict_interval(
                interval=resized_interval,
                requested_outputs={dna_client.OutputType.CONTACT_MAPS},
                organism=organism,
                ontology_terms=None
            )
            wt_contact_map = wt_output.contact_maps
            print(f"   New contact map shape: {wt_contact_map.values.shape}")
            if wt_contact_map.values.shape[-1] > 0:
                print(f"   ✓ Found {wt_contact_map.values.shape[-1]} cell types with contact data")
                print(f"   Using: {wt_contact_map.ontology_terms[0].ontology_curie}")
            else:
                raise ValueError("No contact map data available for this region")

    except Exception as e:
        print(f"   ✗ Error predicting wild-type: {e}")
        return

    print(f"\n2. Predicting DELETION VARIANT contact map...")
    try:
        # Note: AlphaGenome's variant API may not support large structural variants
        # This is an experimental feature - if it fails, we'll note the limitation
        del_output = dna_model.predict_variant(
            interval=resized_interval,
            variant=deletion_variant,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            organism=organism,
            ontology_terms=None  # Use same cell types as wild-type
        )
        print("   ✓ Deletion prediction complete")

        # Extract alternate (deletion) contact map
        del_contact_map = del_output.alternate.contact_maps
        print(f"   Deletion contact map shape: {del_contact_map.values.shape}")

        has_deletion_prediction = True

    except Exception as e:
        print(f"   ✗ Error predicting deletion variant: {e}")
        print(f"   Note: Large structural variants may not be fully supported")
        print(f"   Showing wild-type contact map only")
        has_deletion_prediction = False

    # Visualization
    print(f"\n3. Creating visualizations...")

    if has_deletion_prediction:
        # Side-by-side comparison
        fig = plt.figure(figsize=(16, 7))
        gs = GridSpec(1, 3, width_ratios=[1, 1, 0.05], wspace=0.3)

        # Wild-type contact map
        ax1 = fig.add_subplot(gs[0])
        wt_matrix = wt_contact_map.values[0]  # First cell type
        im1 = ax1.imshow(wt_matrix, cmap='Reds', aspect='auto',
                        interpolation='nearest', vmin=0, vmax=1)
        ax1.set_title(f'Wild-Type Contact Map\n{interval_str}', fontsize=12, fontweight='bold')
        ax1.set_xlabel('Genomic Position (bins)')
        ax1.set_ylabel('Genomic Position (bins)')

        # Mark deletion region
        del_start_bin = int((deletion_start - resized_interval.start) / wt_contact_map.resolution)
        del_end_bin = int((deletion_end - resized_interval.start) / wt_contact_map.resolution)
        ax1.axhline(del_start_bin, color='blue', linestyle='--', linewidth=1, alpha=0.5)
        ax1.axhline(del_end_bin, color='blue', linestyle='--', linewidth=1, alpha=0.5)
        ax1.axvline(del_start_bin, color='blue', linestyle='--', linewidth=1, alpha=0.5)
        ax1.axvline(del_end_bin, color='blue', linestyle='--', linewidth=1, alpha=0.5)

        # Deletion contact map
        ax2 = fig.add_subplot(gs[1])
        del_matrix = del_contact_map.values[0]
        im2 = ax2.imshow(del_matrix, cmap='Reds', aspect='auto',
                        interpolation='nearest', vmin=0, vmax=1)
        ax2.set_title(f'After Deletion\n(Removed {deletion_end - deletion_start:,} bp)',
                     fontsize=12, fontweight='bold')
        ax2.set_xlabel('Genomic Position (bins)')
        ax2.set_ylabel('Genomic Position (bins)')

        # Shared colorbar
        cax = fig.add_subplot(gs[2])
        plt.colorbar(im2, cax=cax, label='Contact Probability')

        plt.suptitle(f'{title_prefix} Effect on Chromatin Contact Maps',
                    fontsize=14, fontweight='bold', y=0.98)

    else:
        # Wild-type only
        fig, ax = plt.subplots(figsize=(8, 7))
        wt_matrix = wt_contact_map.values[0]
        im = ax.imshow(wt_matrix, cmap='Reds', aspect='auto',
                      interpolation='nearest', vmin=0, vmax=1)
        ax.set_title(f'Wild-Type Contact Map\n{interval_str}', fontsize=12, fontweight='bold')
        ax.set_xlabel('Genomic Position (bins)')
        ax.set_ylabel('Genomic Position (bins)')

        # Mark deletion region
        del_start_bin = int((deletion_start - resized_interval.start) / wt_contact_map.resolution)
        del_end_bin = int((deletion_end - resized_interval.start) / wt_contact_map.resolution)
        ax.axhline(del_start_bin, color='blue', linestyle='--', linewidth=1.5, alpha=0.7,
                  label=f'Deletion region ({deletion_end - deletion_start:,} bp)')
        ax.axhline(del_end_bin, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(del_start_bin, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.axvline(del_end_bin, color='blue', linestyle='--', linewidth=1.5, alpha=0.7)
        ax.legend(loc='upper right')

        plt.colorbar(im, ax=ax, label='Contact Probability')

    # Save figure
    safe_interval = interval_str.replace(':', '_').replace('-', '_')
    output_path = f"media/deletion_effect_{safe_interval}_{deletion_start}_{deletion_end}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"   ✓ Saved to {output_path}")
    plt.close()

    print(f"\n{'='*70}")
    print("Analysis complete!")
    print(f"{'='*70}\n")


def main():
    """
    Run deletion effect analysis on mouse genomic regions.

    We'll test deletions of different genomic features:
    1. TAD boundary deletion (predicted to disrupt TAD structure)
    2. CTCF binding site deletion (predicted to weaken loop anchoring)
    3. Regulatory element deletion

    AlphaGenome supports both human and mouse genomes!
    """

    # Example 1: Mouse Sox11 region (chr12) - TAD boundary region
    print("\n" + "="*70)
    print("EXPERIMENT 1: TAD Boundary Deletion in Sox11 Region (Mouse)")
    print("="*70)

    sox11_interval = "chr12:26000000-27000000"
    # Delete ~50kb in the middle (hypothetical TAD boundary)
    sox11_del_start = 26500000
    sox11_del_end = 26550000

    visualize_deletion_effect(
        interval_str=sox11_interval,
        deletion_start=sox11_del_start,
        deletion_end=sox11_del_end,
        cell_type=None,  # Auto-detect available cell types
        organism=dna_client.Organism.MUS_MUSCULUS,
        title_prefix="Sox11 TAD Boundary Deletion"
    )

    # Example 2: Mouse Mir9-2 region (chr13) - regulatory region
    print("\n" + "="*70)
    print("EXPERIMENT 2: Regulatory Element Deletion in Mir9-2 Region (Mouse)")
    print("="*70)

    mir9_interval = "chr13:83500000-84500000"
    # Delete a smaller regulatory region
    mir9_del_start = 84000000
    mir9_del_end = 84020000

    visualize_deletion_effect(
        interval_str=mir9_interval,
        deletion_start=mir9_del_start,
        deletion_end=mir9_del_end,
        cell_type=None,  # Auto-detect available cell types
        organism=dna_client.Organism.MUS_MUSCULUS,
        title_prefix="Mir9-2 Regulatory Element Deletion"
    )

    print("\n" + "="*70)
    print("All deletion analyses complete!")
    print("Check the 'media/' directory for visualizations")
    print("="*70)


if __name__ == "__main__":
    main()

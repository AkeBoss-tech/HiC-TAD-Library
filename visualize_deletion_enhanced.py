#!/usr/bin/env python3
"""
Enhanced DNA deletion effect visualization with gene annotations and TAD structures.

Creates publication-quality contact map visualizations similar to AlphaGenome papers,
showing gene annotations, contact maps, and deletion effects side-by-side.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from dotenv import load_dotenv

from alphagenome.data import genome, gene_annotation, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components
import pandas as pd

# Load environment variables
load_dotenv()
API_KEY = os.getenv('ALPHA_GENOME_API_KEY')

# Initialize AlphaGenome model
print("Initializing AlphaGenome model...")
dna_model = dna_client.create(API_KEY)

# Load gene annotations
print("Loading gene annotations...")
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf = pd.read_feather(gtf_url)

# Filter to protein-coding genes
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1', '2']
)

# Create transcript extractor
transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
print("Caching transcripts...")
transcript_extractor.cache_transcripts()

# Create longest transcript extractor
gtf_longest = gene_annotation.filter_to_longest_transcript(gtf_transcript)
longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest)
longest_transcript_extractor.cache_transcripts()

print("✓ Ready!")


def visualize_enhanced_deletion(
    interval_str,
    deletion_start,
    deletion_end,
    cell_type='EFO:0002824',  # HCT116 colon carcinoma
    organism=dna_client.Organism.HOMO_SAPIENS,
    title="Contact Map Analysis"
):
    """
    Create enhanced publication-quality visualization of deletion effects.

    Features:
    - Gene annotations at the top
    - Wild-type and deletion contact maps side-by-side
    - Proper genomic coordinates
    - TAD structure visualization
    - Deletion region highlighted
    """
    print(f"\n{'='*80}")
    print(f"Enhanced Deletion Analysis: {interval_str}")
    print(f"Deletion: {deletion_start:,} - {deletion_end:,} ({deletion_end - deletion_start:,} bp)")
    print(f"{'='*80}\n")

    # Parse and resize interval
    interval = genome.Interval.from_str(interval_str)
    resized_interval = interval.resize(2**20)  # 1,048,576 bp
    print(f"Analysis interval: {resized_interval}")

    # Create deletion variant
    deletion_variant = genome.Variant(
        chromosome=interval.chromosome,
        position=deletion_start - 1,
        reference_bases='N' * (deletion_end - deletion_start + 1),
        alternate_bases='N',
    )

    # Get gene annotations
    print("\n1. Fetching gene annotations...")
    transcripts = longest_transcript_extractor.extract(resized_interval)
    print(f"   Found {len(transcripts)} transcripts")

    # Predict wild-type
    print("\n2. Predicting wild-type contact map...")
    wt_output = dna_model.predict_interval(
        interval=resized_interval,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=[cell_type],
        organism=organism
    )
    print(f"   ✓ Contact map shape: {wt_output.contact_maps.values.shape}")

    # Predict deletion variant
    print("\n3. Predicting deletion variant contact map...")
    try:
        del_output = dna_model.predict_variant(
            interval=resized_interval,
            variant=deletion_variant,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            ontology_terms=[cell_type],
            organism=organism
        )
        print(f"   ✓ Deletion contact map shape: {del_output.alternate.contact_maps.values.shape}")
        has_deletion = True
    except Exception as e:
        print(f"   ✗ Deletion prediction failed: {e}")
        has_deletion = False

    # Create enhanced visualization
    print("\n4. Creating publication-quality visualization...")

    if has_deletion:
        # Use AlphaGenome's plotting library for professional appearance
        fig = plt.figure(figsize=(18, 10))

        # Top: Gene annotations spanning both contact maps
        ax_genes = plt.subplot2grid((2, 2), (0, 0), colspan=2, fig=fig)

        # Plot gene annotations
        for tx in transcripts:
            y_pos = 0.5
            tx_int = tx.transcript_interval
            # Draw gene body
            ax_genes.plot([tx_int.start, tx_int.end], [y_pos, y_pos],
                         color='darkblue', linewidth=2)
            # Add gene name
            gene_center = (tx_int.start + tx_int.end) / 2
            gene_name = tx.info.get('gene_name', '')
            if gene_name:
                ax_genes.text(gene_center, y_pos + 0.1, gene_name,
                            ha='center', va='bottom', fontsize=8,
                            fontweight='bold')

        # Mark deletion region
        ax_genes.axvspan(deletion_start, deletion_end, alpha=0.3,
                        color='red', label='Deletion')
        ax_genes.set_xlim(resized_interval.start, resized_interval.end)
        ax_genes.set_ylim(0, 1)
        ax_genes.set_ylabel('Genes', fontweight='bold')
        ax_genes.set_title(title, fontsize=14, fontweight='bold', pad=20)
        ax_genes.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
        ax_genes.spines['top'].set_visible(False)
        ax_genes.spines['right'].set_visible(False)
        ax_genes.spines['left'].set_visible(False)
        ax_genes.set_yticks([])

        # Format x-axis
        ax_genes.ticklabel_format(style='plain', axis='x')
        ax_genes.set_xlabel('')

        # Bottom left: Wild-type contact map
        ax_wt = plt.subplot2grid((2, 2), (1, 0), fig=fig)
        wt_matrix = wt_output.contact_maps.values[:, :, 0]  # First cell type

        extent = [resized_interval.start, resized_interval.end,
                 resized_interval.end, resized_interval.start]

        im_wt = ax_wt.imshow(wt_matrix, cmap='Reds', aspect='auto',
                            interpolation='nearest', extent=extent,
                            vmin=0, vmax=np.percentile(wt_matrix, 99))

        # Mark deletion region
        ax_wt.axhline(deletion_start, color='blue', linestyle='--',
                     linewidth=1.5, alpha=0.7)
        ax_wt.axhline(deletion_end, color='blue', linestyle='--',
                     linewidth=1.5, alpha=0.7)
        ax_wt.axvline(deletion_start, color='blue', linestyle='--',
                     linewidth=1.5, alpha=0.7)
        ax_wt.axvline(deletion_end, color='blue', linestyle='--',
                     linewidth=1.5, alpha=0.7)

        ax_wt.set_title('Wild-Type', fontweight='bold', fontsize=12)
        ax_wt.set_ylabel('Chromosome position', fontweight='bold')
        ax_wt.set_xlabel('Chromosome position', fontweight='bold')
        ax_wt.ticklabel_format(style='plain')

        # Add cell type label on left
        cell_label = f"{cell_type}\nWild-Type"
        ax_wt.text(-0.15, 0.5, cell_label, transform=ax_wt.transAxes,
                  rotation=90, va='center', ha='center', fontsize=10,
                  fontweight='bold')

        # Bottom right: Deletion contact map
        ax_del = plt.subplot2grid((2, 2), (1, 1), fig=fig)
        del_matrix = del_output.alternate.contact_maps.values[:, :, 0]

        im_del = ax_del.imshow(del_matrix, cmap='Reds', aspect='auto',
                              interpolation='nearest', extent=extent,
                              vmin=0, vmax=np.percentile(del_matrix, 99))

        ax_del.set_title(f'After Deletion ({deletion_end - deletion_start:,} bp)',
                        fontweight='bold', fontsize=12)
        ax_del.set_ylabel('Chromosome position', fontweight='bold')
        ax_del.set_xlabel('Chromosome position', fontweight='bold')
        ax_del.ticklabel_format(style='plain')

        # Add cell type label on left
        cell_label_del = f"{cell_type}\nDeletion"
        ax_del.text(-0.15, 0.5, cell_label_del, transform=ax_del.transAxes,
                   rotation=90, va='center', ha='center', fontsize=10,
                   fontweight='bold')

        # Add shared colorbar
        cbar = plt.colorbar(im_del, ax=[ax_wt, ax_del],
                           orientation='vertical', pad=0.02, aspect=30)
        cbar.set_label('Contact Probability', fontweight='bold')

        plt.tight_layout()

    else:
        # Wild-type only if deletion failed
        fig, (ax_genes, ax_wt) = plt.subplots(2, 1, figsize=(14, 10),
                                              gridspec_kw={'height_ratios': [1, 4]})

        # Gene annotations
        for tx in transcripts:
            y_pos = 0.5
            tx_int = tx.transcript_interval
            ax_genes.plot([tx_int.start, tx_int.end], [y_pos, y_pos],
                         color='darkblue', linewidth=2)
            gene_center = (tx_int.start + tx_int.end) / 2
            gene_name = tx.info.get('gene_name', '')
            if gene_name:
                ax_genes.text(gene_center, y_pos + 0.1, gene_name,
                            ha='center', va='bottom', fontsize=8,
                            fontweight='bold')

        ax_genes.axvspan(deletion_start, deletion_end, alpha=0.3,
                        color='red', label='Deletion')
        ax_genes.set_xlim(resized_interval.start, resized_interval.end)
        ax_genes.set_ylim(0, 1)
        ax_genes.set_ylabel('Genes', fontweight='bold')
        ax_genes.set_title(title, fontsize=14, fontweight='bold')
        ax_genes.legend(loc='upper right', frameon=True, fancybox=True, shadow=True)
        ax_genes.spines['top'].set_visible(False)
        ax_genes.spines['right'].set_visible(False)
        ax_genes.spines['left'].set_visible(False)
        ax_genes.set_yticks([])

        # Contact map
        wt_matrix = wt_output.contact_maps.values[:, :, 0]
        extent = [resized_interval.start, resized_interval.end,
                 resized_interval.end, resized_interval.start]

        im = ax_wt.imshow(wt_matrix, cmap='Reds', aspect='auto',
                         interpolation='nearest', extent=extent,
                         vmin=0, vmax=np.percentile(wt_matrix, 99))

        ax_wt.set_title('Wild-Type Contact Map', fontweight='bold')
        ax_wt.set_ylabel('Chromosome position', fontweight='bold')
        ax_wt.set_xlabel('Chromosome position', fontweight='bold')
        ax_wt.ticklabel_format(style='plain')

        plt.colorbar(im, ax=ax_wt, label='Contact Probability')
        plt.tight_layout()

    # Save figure
    safe_interval = interval_str.replace(':', '_').replace('-', '_')
    output_path = f"media/enhanced_deletion_{safe_interval}_{deletion_start}_{deletion_end}.png"
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"\n   ✓ Saved to {output_path}")
    plt.close()

    print(f"\n{'='*80}")
    print("Enhanced visualization complete!")
    print(f"{'='*80}\n")


def main():
    """
    Create enhanced deletion visualizations for key genomic regions.
    """

    print("\n" + "="*80)
    print("ENHANCED DELETION ANALYSIS WITH GENE ANNOTATIONS")
    print("="*80)

    # Example 1: OCT4 region - pluripotency master regulator
    print("\n📍 Region 1: OCT4 (Pluripotency)")
    visualize_enhanced_deletion(
        interval_str="chr6:30630000-31630000",
        deletion_start=31100000,
        deletion_end=31150000,
        cell_type='EFO:0002824',  # HCT116
        organism=dna_client.Organism.HOMO_SAPIENS,
        title="Predicted TADs near OCT4 - TAD Boundary Deletion Effect"
    )

    # Example 2: SOX2 region - neural differentiation
    print("\n📍 Region 2: SOX2 (Neural Development)")
    visualize_enhanced_deletion(
        interval_str="chr3:181210000-182210000",
        deletion_start=181700000,
        deletion_end=181750000,
        cell_type='EFO:0002824',  # HCT116
        organism=dna_client.Organism.HOMO_SAPIENS,
        title="Predicted TADs near SOX2 - Regulatory Element Deletion"
    )

    # Example 3: NANOG region - pluripotency
    print("\n📍 Region 3: NANOG (Pluripotency)")
    visualize_enhanced_deletion(
        interval_str="chr12:7280000-8280000",
        deletion_start=7750000,
        deletion_end=7800000,
        cell_type='EFO:0002824',  # HCT116
        organism=dna_client.Organism.HOMO_SAPIENS,
        title="Predicted TADs near NANOG - Enhancer Deletion"
    )

    print("\n" + "="*80)
    print("✅ All enhanced visualizations complete!")
    print("Check the 'media/' directory for high-resolution figures")
    print("="*80)


if __name__ == "__main__":
    main()

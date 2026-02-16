import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

from alphagenome.data import genome, gene_annotation, track_data, transcript
from alphagenome.models import dna_client
from alphagenome.visualization import plot_components

# 1. Setup and Environment
load_dotenv()
api_key = os.getenv("ALPHA_GENOME_API_KEY")

if not api_key:
    raise ValueError("ALPHA_GENOME_API_KEY not found in .env file")

# Initialize DNA client
dna_model = dna_client.create(api_key)

# 2. Reference Data Setup
print("Loading reference data...")
# Load gene annotations (from GENCODE hg38)
# Using the feather format recommended in the documentation for speed
gtf_url = 'https://storage.googleapis.com/alphagenome/reference/gencode/hg38/gencode.v46.annotation.gtf.gz.feather'
gtf = pd.read_feather(gtf_url)

# Filter to protein-coding genes and highly supported transcripts (level 1)
gtf_transcript = gene_annotation.filter_transcript_support_level(
    gene_annotation.filter_protein_coding(gtf), ['1']
)

# Extractor for identifying transcripts in a region
transcript_extractor = transcript.TranscriptExtractor(gtf_transcript)
print("Caching transcripts (this may take a minute)...")
transcript_extractor.cache_transcripts()

# Extractor for only the longest transcript per gene (cleaner for some plots)
gtf_longest_transcript = gene_annotation.filter_to_longest_transcript(gtf_transcript)
longest_transcript_extractor = transcript.TranscriptExtractor(gtf_longest_transcript)
longest_transcript_extractor.cache_transcripts()

# 3. Define Regions of Interest
# APOL4 region on chr22 (hg38)
apol4_interval = 'chr22:36,245,000-36,275,000'

# 4. Visualization Functions

def visualize_gene_expression(interval_str, title="Gene Expression and Tracks"):
    print(f"Visualizing Gene Expression for {interval_str}...")
    interval = genome.Interval.from_str(interval_str)
    transcripts = transcript_extractor.extract(interval)

    # Resize interval to model's supported length (2^20 = 1,048,576 bp)
    resized_interval = interval.resize(2**20)

    # Use a common cell type ontology term (K562 cell line - chronic myelogenous leukemia)
    ontology_terms = ['EFO:0002067']

    # Predict across specified cell types
    output = dna_model.predict_interval(
        interval=resized_interval,
        requested_outputs={dna_client.OutputType.RNA_SEQ},
        ontology_terms=ontology_terms
    )

    # Plotting
    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(
                tdata=output.rna_seq,
                ylabel_template='{biosample_name}',
                filled=True,
            ),
        ],
        interval=resized_interval,
        title=title,
        despine_keep_bottom=True,
    )

    output_path = f"media/alphagenome_expression_{interval_str.replace(':', '_').replace('-', '_')}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()

def visualize_contact_maps(interval_str, title="Predicted Contact Maps (TADs)"):
    print(f"Visualizing Contact Maps for {interval_str}...")
    interval = genome.Interval.from_str(interval_str)

    # Resize interval to model's supported length (2^20 = 1,048,576 bp)
    resized_interval = interval.resize(2**20)

    # Contact maps are usually predicted for specific cell lines
    # HCT116 colon carcinoma cell line (EFO:0002824)
    ontology_terms = ['EFO:0002824']

    output = dna_model.predict_interval(
        interval=resized_interval,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=ontology_terms,
    )
    
    longest_transcripts = longest_transcript_extractor.extract(resized_interval)

    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.ContactMaps(
                tdata=output.contact_maps,
                ylabel_template='{biosample_name}\n{name}',
                cmap='autumn_r',
                vmax=1.0,
            ),
        ],
        interval=resized_interval,
        title=title,
    )

    output_path = f"media/alphagenome_contacts_{interval_str.replace(':', '_').replace('-', '_')}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()

def visualize_virtual_4c(interval_str, anchor_pos, title="Virtual 4C"):
    print(f"Visualizing Virtual 4C for {interval_str} with anchor at {anchor_pos}...")
    interval = genome.Interval.from_str(interval_str)

    # Resize interval to model's supported length (2^20 = 1,048,576 bp)
    resized_interval = interval.resize(2**20)

    ontology_terms = ['EFO:0002824']
    output = dna_model.predict_interval(
        interval=resized_interval,
        requested_outputs={dna_client.OutputType.CONTACT_MAPS},
        ontology_terms=ontology_terms,
    )
    
    # Convert contact map to virtual 4C track from the given anchor
    v4c_track = output.contact_maps.to_virtual_4c(anchor_pos)

    longest_transcripts = longest_transcript_extractor.extract(resized_interval)

    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(longest_transcripts),
            plot_components.Tracks(
                tdata=v4c_track,
                ylabel_template='Virtual 4C\n(Anchor: {anchor})',
                filled=True,
                color='darkblue'
            ),
        ],
        interval=resized_interval,
        title=title,
        annotations=[plot_components.IntervalAnnotation([f"{resized_interval.chromosome}:{anchor_pos}-{anchor_pos+2048}"], color='red', alpha=0.5)]
    )

    output_path = f"media/alphagenome_v4c_{interval_str.replace(':', '_').replace('-', '_')}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()

def visualize_ctcf_signal(interval_str, title="Predicted CTCF Binding"):
    print(f"Visualizing CTCF signal for {interval_str}...")
    interval = genome.Interval.from_str(interval_str)

    # Resize interval to model's supported length (2^20 = 1,048,576 bp)
    resized_interval = interval.resize(2**20)

    # Use K562 cell line for CTCF predictions
    ontology_terms = ['EFO:0002067']

    # Fetch ChIP-TF predictions
    output = dna_model.predict_interval(
        interval=resized_interval,
        requested_outputs={dna_client.OutputType.CHIP_TF},
        ontology_terms=ontology_terms
    )
    
    # Filter for CTCF and average across tissues
    ctcf_mask = output.chip_tf.metadata['transcription_factor'] == 'CTCF'
    mean_ctcf = output.chip_tf.values[:, ctcf_mask].mean(axis=1)
    
    # Construct a new TrackData object for the mean
    tdata_mean_ctcf = track_data.TrackData(
        values=mean_ctcf[:, None],
        metadata=pd.DataFrame({
            'transcription_factor': ['CTCF'],
            'name': ['mean'],
            'strand': ['.']
        }),
        interval=output.chip_tf.interval,
        resolution=output.chip_tf.resolution,
    )
    
    transcripts = transcript_extractor.extract(resized_interval)

    plot_components.plot(
        [
            plot_components.TranscriptAnnotation(transcripts),
            plot_components.Tracks(
                tdata=tdata_mean_ctcf,
                ylabel_template='{name} {transcription_factor}',
                filled=True
            ),
        ],
        interval=resized_interval,
        title=title,
        despine_keep_bottom=True,
    )

    output_path = f"media/alphagenome_ctcf_{interval_str.replace(':', '_').replace('-', '_')}.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()

if __name__ == "__main__":
    if not os.path.exists("media"):
        os.makedirs("media")
        
    print("Starting AlphaGenome Visualization Modality Tour...")
    
    # 5. Defined Regions for Exploration
    regions = {
        # Human genome (hg38) regions
        "APOL4": "chr22:36245000-36275000",
        "OCT4": "chr6:30630000-31630000", # Pluripotency
        "NANOG": "chr12:7280000-8280000", # Pluripotency
        "SOX2": "chr3:181210000-182210000", # Pluripotency
        "PAX6": "chr11:31300000-32300000", # Neural Differentiation
        "NKX2-5": "chr5:173100000-174100000", # Cardiac Differentiation

        # Mouse genome regions (may need conversion to hg38 if AlphaGenome only supports human)
        # Note: AlphaGenome is trained on human genome - these may not work
        "Sox11_Mouse": "chr12:26000000-28000000", # Mouse Sox11 region
        "Mir9-2_Mouse": "chr13:83500000-84500000" # Mouse Mir9-2 region
    }

    # Define which regions to analyze (you can add/remove regions from this list)
    # Options: "APOL4", "OCT4", "NANOG", "SOX2", "PAX6", "NKX2-5", "Sox11_Mouse", "Mir9-2_Mouse"
    regions_to_analyze = ["Sox11_Mouse", "Mir9-2_Mouse"]

    # Run visualizations for each selected region
    for target_gene in regions_to_analyze:
        selected_interval = regions[target_gene]

        print(f"\n--- Analyzing {target_gene} region: {selected_interval} ---")

        try:
            # Gene expression for the specific gene window
            visualize_gene_expression(selected_interval, f"{target_gene} Gene Expression")

            # For contact maps and v4c, we usually want a slightly larger window (~500kb-1Mb)
            # The regions defined above are already ~1Mb chunks
            visualize_contact_maps(selected_interval, f"Predicted TADs near {target_gene}")

            # Virtual 4C with anchor at the center (approximate promoter location)
            # NOTE: to_virtual_4c method not available in current AlphaGenome API
            # interval_obj = genome.Interval.from_str(selected_interval)
            # anchor = (interval_obj.start + interval_obj.end) // 2
            # visualize_virtual_4c(selected_interval, anchor, f"Virtual 4C for {target_gene}")

            # CTCF Binding
            visualize_ctcf_signal(selected_interval, f"Mean Predicted CTCF Binding at {target_gene}")

            print(f"✓ Successfully completed visualizations for {target_gene}")
        except Exception as e:
            print(f"✗ Error analyzing {target_gene}: {str(e)}")
            print(f"  This may occur if the region is not in the supported genome assembly.")
            continue

    print("\nVisualization Tour Complete! Check the 'media/' directory for results.")

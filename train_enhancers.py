"""
Train the EnhancerCNN using AlphaGenome ATAC-seq predictions as proxy labels.

Strategy
--------
- Make one AlphaGenome call per 1 Mb window (3 calls total for 2 regions)
- Extract predicted motor-neuron ATAC signal at each candidate bin's position
- Train the CNN on (500 bp sequence, ATAC score) pairs
- Save model to data/processed/enhancer_cnn.pt
- Re-run visualize_enhancers.py to produce CNN-scored figures

Cell type: CL:0000100 (motor neuron ATAC) — best available proxy for
           Sox11 (neural TF) and Mir9-2 (neural miRNA) regulatory context.
Organism:  MUS_MUSCULUS (mm10)
"""

import os
import sys
import cooler
import numpy as np
import pandas as pd
from dotenv import load_dotenv

from src.tad_boundaries import (
    score_boundary_prominence, call_tad_intervals, multiscale_insulation
)
from src.enhancer_id import (
    extract_candidate_enhancers, fetch_sequences, train_enhancer_cnn
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

FILE_PATH  = "data/raw/mouse_microc.mcool"
FASTA_PATH = "data/raw/mm10.fa"
MODEL_PATH = "data/processed/enhancer_cnn.pt"
RESOLUTION = 5000

# AlphaGenome cell type — motor neuron ATAC is the closest neural proxy available
ATAC_ONTOLOGY = "CL:0000100"

regions = {
    "Sox11_Chr12_a": "chr12:26,000,000-27,000,000",
    "Sox11_Chr12_b": "chr12:27,000,000-28,000,000",
    "Mir9-2_Chr13":  "chr13:83,500,000-84,500,000",
}


# ---------------------------------------------------------------------------
# AlphaGenome ATAC extraction
# ---------------------------------------------------------------------------

def get_atac_for_region(model, coordinates: str) -> tuple:
    """
    Query AlphaGenome for ATAC predictions across a ~1 Mb region.

    Returns (positions_bp, atac_signal) where positions_bp is a numpy
    array of genomic positions (bp) at 128 bp resolution.
    """
    from alphagenome.data import genome
    from alphagenome.models import dna_model, dna_output

    chrom, rest = coordinates.split(":")
    start, end = [int(x.replace(",", "")) for x in rest.split("-")]
    center = (start + end) // 2

    # AlphaGenome requires exactly 2^20 bp window
    interval = genome.Interval(chrom, center, center).resize(2**20)

    print(f"  AlphaGenome ATAC query: {chrom}:{interval.start}-{interval.end}")
    output = model.predict_interval(
        interval,
        organism=dna_model.Organism.MUS_MUSCULUS,
        requested_outputs=[dna_output.OutputType.ATAC],
        ontology_terms=[ATAC_ONTOLOGY],
    )

    # output.atac is shape (n_positions, 1) at 128 bp resolution
    atac_values = output.atac.values.squeeze()   # (n_positions,)
    n = len(atac_values)
    positions = np.array([interval.start + i * 128 + 64 for i in range(n)])

    return positions, atac_values


def assign_atac_to_candidates(
    candidates: pd.DataFrame,
    atac_positions: np.ndarray,
    atac_values: np.ndarray,
) -> np.ndarray:
    """
    For each candidate bin, find the nearest ATAC prediction position
    and return its signal value.
    """
    labels = np.zeros(len(candidates), dtype=np.float32)
    bin_mids = ((candidates["start"] + candidates["end"]) / 2).values

    for i, mid in enumerate(bin_mids):
        idx = np.argmin(np.abs(atac_positions - mid))
        labels[i] = float(atac_values[idx])

    # Normalize to [0, 1]
    lo, hi = labels.min(), labels.max()
    if hi > lo:
        labels = (labels - lo) / (hi - lo)
    return labels


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    load_dotenv()
    api_key = os.getenv("ALPHA_GENOME_API_KEY")
    if not api_key:
        print("ERROR: ALPHA_GENOME_API_KEY not found in .env")
        sys.exit(1)

    print("=== Enhancer CNN Training ===")
    print(f"ATAC proxy cell type: {ATAC_ONTOLOGY} (motor neuron)")
    print(f"Organism: MUS_MUSCULUS (mm10)")

    # Load AlphaGenome
    from alphagenome.models import dna_client
    ag_model = dna_client.create(api_key)

    # Load cooler
    uri = f"{FILE_PATH}::resolutions/{RESOLUTION}"
    clr = cooler.Cooler(uri)

    all_sequences = []
    all_labels    = []

    for region_name, coordinates in regions.items():
        chrom, rest = coordinates.split(":")
        start, end = [int(x.replace(",", "")) for x in rest.split("-")]
        print(f"\n--- {region_name} ({coordinates}) ---")

        # TAD pipeline
        ins_table, windows = multiscale_insulation(
            clr, coordinates, window_sizes_bp=[50_000]
        )
        boundary_df = score_boundary_prominence(ins_table, 50_000)
        tad_df = call_tad_intervals(boundary_df, clr)
        print(f"  TADs called: {len(tad_df)}")

        # Candidate extraction
        candidates = extract_candidate_enhancers(clr, coordinates, tad_df, boundary_df)
        print(f"  Candidates: {len(candidates)}")

        if len(candidates) == 0:
            print("  Skipping (no candidates).")
            continue

        # AlphaGenome ATAC labels
        print("  Querying AlphaGenome ATAC...")
        try:
            atac_pos, atac_vals = get_atac_for_region(ag_model, coordinates)
            labels = assign_atac_to_candidates(candidates, atac_pos, atac_vals)
            print(f"  ATAC labels: min={labels.min():.3f}  max={labels.max():.3f}  mean={labels.mean():.3f}")
        except Exception as e:
            print(f"  AlphaGenome failed ({e}), using contact score as fallback label.")
            raw = candidates["contact_score"].values.astype(np.float32)
            lo, hi = raw.min(), raw.max()
            labels = (raw - lo) / (hi - lo) if hi > lo else raw

        # Sequence extraction
        print("  Fetching sequences from mm10.fa...")
        sequences = fetch_sequences(candidates, FASTA_PATH, window_bp=500)
        print(f"  Sequences shape: {sequences.shape}")

        all_sequences.append(sequences)
        all_labels.append(labels)

    if not all_sequences:
        print("No training data collected. Exiting.")
        sys.exit(1)

    X = np.concatenate(all_sequences, axis=0)
    y = np.concatenate(all_labels,    axis=0)
    print(f"\n=== Training CNN on {len(X)} examples ===")
    print(f"Label distribution: min={y.min():.3f}  max={y.max():.3f}  mean={y.mean():.3f}")

    model = train_enhancer_cnn(
        X, y,
        epochs=30,
        lr=1e-3,
        batch_size=32,
        save_path=MODEL_PATH,
        verbose=True,
    )
    print(f"\nModel saved to {MODEL_PATH}")

    # Re-run visualization with trained model
    print("\n=== Re-running visualize_enhancers.py with trained model ===")
    import subprocess, sys as _sys
    result = subprocess.run(
        [_sys.executable, "visualize_enhancers.py"],
        capture_output=False
    )
    if result.returncode == 0:
        print("\nDone. Check media/*_enhancers.png for CNN-scored figures.")
    else:
        print("\nvisualize_enhancers.py exited with an error — check output above.")


if __name__ == "__main__":
    main()

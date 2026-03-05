"""
Synthetic enhancer design script.

Uses a trained EnhancerCNN to generate novel sequences predicted to have
high enhancer activity in the Sox11 / Mir9-2 neural regulatory context.

Two strategies are run:

1. Gradient-based optimisation — 10 sequences designed from random start.
2. Motif-insertion hill-climb — each of the top-5 real candidate bins is
   refined by greedily inserting neural TF motifs (Sox11, AP-1, NeuroD1, etc.).

Results are saved to:
  data/processed/gradient_designs.tsv   — gradient-optimised sequences
  data/processed/motif_designs.tsv      — motif-insertion designs

Usage
-----
    python design_enhancers.py

Requires
--------
- data/processed/enhancer_cnn.pt  (trained model from train_enhancers.py)
- data/raw/mm10.fa                 (for fetching real seed sequences)
- data/raw/mouse_microc.mcool
"""

import os
import cooler
import numpy as np
import pandas as pd

from src.tad_boundaries import (
    score_boundary_prominence, call_tad_intervals, multiscale_insulation
)
from src.enhancer_id import (
    extract_candidate_enhancers, fetch_sequences, EnhancerCNN
)
from src.enhancer_design import (
    gradient_design, motif_insert_design, score_sequences, summarise_designs,
    NEURAL_MOTIFS, gc_content,
)

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------

FILE_PATH  = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
FASTA_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mm10.fa"
MODEL_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/processed/enhancer_cnn.pt"
OUT_DIR    = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/processed"
RESOLUTION = 5000

REGIONS = {
    "Sox11_Chr12": "chr12:26,000,000-28,000,000",
    "Mir9-2_Chr13": "chr13:83,500,000-84,500,000",
}

N_GRADIENT_DESIGNS = 10   # per region
N_MOTIF_SEEDS      = 5    # top-N real candidates to refine
N_MOTIF_ROUNDS     = 3    # hill-climb rounds


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _run_tad_pipeline(clr, coordinates):
    ins_table, _ = multiscale_insulation(clr, coordinates, window_sizes_bp=[50_000])
    boundary_df = score_boundary_prominence(ins_table, 50_000)
    tad_df = call_tad_intervals(boundary_df, clr)
    return tad_df, boundary_df


def _get_candidates_with_scores(clr, coordinates, model):
    """Return candidates sorted by CNN score (descending)."""
    tad_df, boundary_df = _run_tad_pipeline(clr, coordinates)
    candidates = extract_candidate_enhancers(clr, coordinates, tad_df, boundary_df)
    if len(candidates) == 0:
        return candidates

    sequences = fetch_sequences(candidates, FASTA_PATH, window_bp=500)
    cnn_scores = score_sequences(model, [
        _tensor_to_str(sequences[i]) for i in range(len(sequences))
    ])
    candidates = candidates.copy()
    candidates["cnn_score"] = cnn_scores
    candidates = candidates.sort_values("cnn_score", ascending=False).reset_index(drop=True)
    return candidates, sequences


def _tensor_to_str(arr: np.ndarray) -> str:
    """Convert (4, L) one-hot numpy array → ACGT string."""
    bases = ['A', 'C', 'G', 'T']
    indices = arr.argmax(axis=0)
    return ''.join(bases[i] for i in indices)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    if not os.path.exists(MODEL_PATH):
        print(f"ERROR: Trained model not found at {MODEL_PATH}")
        print("Run train_enhancers.py first.")
        return

    if not os.path.exists(FASTA_PATH):
        print(f"ERROR: Reference genome not found at {FASTA_PATH}")
        print("Download mm10.fa and index with samtools faidx.")
        return

    print("=== Synthetic Enhancer Design ===")
    print(f"Model: {MODEL_PATH}")

    # Load trained CNN
    cnn = EnhancerCNN()
    cnn.load(MODEL_PATH)
    print("CNN loaded.")

    # Load cooler
    uri = f"{FILE_PATH}::resolutions/{RESOLUTION}"
    clr = cooler.Cooler(uri)

    all_gradient_rows = []
    all_motif_rows    = []

    for region_name, coordinates in REGIONS.items():
        print(f"\n{'='*60}")
        print(f"Region: {region_name}  ({coordinates})")
        print(f"{'='*60}")

        # Get real candidates, scored by CNN
        tad_df, boundary_df = _run_tad_pipeline(clr, coordinates)
        candidates = extract_candidate_enhancers(clr, coordinates, tad_df, boundary_df)

        if len(candidates) == 0:
            print("No candidates found — skipping.")
            continue

        sequences_arr = fetch_sequences(candidates, FASTA_PATH, window_bp=500)
        seed_seqs = [_tensor_to_str(sequences_arr[i]) for i in range(len(sequences_arr))]
        cnn_scores = score_sequences(cnn, seed_seqs)
        candidates = candidates.copy()
        candidates["cnn_score"] = cnn_scores
        # Track original position so seed_seqs[orig_idx] stays aligned after sort
        candidates["_orig_idx"] = range(len(candidates))
        candidates = candidates.sort_values("cnn_score", ascending=False).reset_index(drop=True)

        print(f"Candidates: {len(candidates)}  |  top CNN score: {cnn_scores.max():.4f}")

        # ----------------------------------------------------------------
        # Strategy 1: Gradient-based design
        # ----------------------------------------------------------------
        print(f"\n--- Gradient design ({N_GRADIENT_DESIGNS} sequences) ---")
        grad_seqs, grad_scores = gradient_design(
            cnn,
            n_sequences=N_GRADIENT_DESIGNS,
            length=500,
            steps=600,
            lr=0.05,
            gc_penalty=0.5,
            verbose=True,
        )
        for seq, sc in zip(grad_seqs, grad_scores):
            all_gradient_rows.append({
                "region":     region_name,
                "strategy":   "gradient",
                "cnn_score":  round(sc, 4),
                "gc_content": round(gc_content(seq), 4),
                "length":     len(seq),
                "sequence":   seq,
            })

        # ----------------------------------------------------------------
        # Strategy 2: Motif-insertion from top real candidates
        # ----------------------------------------------------------------
        print(f"\n--- Motif-insertion design (top {N_MOTIF_SEEDS} seeds) ---")
        n_seeds = min(N_MOTIF_SEEDS, len(candidates))
        for idx in range(n_seeds):
            row = candidates.iloc[idx]
            orig_idx = int(row["_orig_idx"])
            seed_seq = seed_seqs[orig_idx]
            seed_score = float(row["cnn_score"])
            print(f"\n  Seed {idx+1}: {row['chrom']}:{row['start']}-{row['end']}  "
                  f"score={seed_score:.4f}")
            best_seq, best_score, history = motif_insert_design(
                cnn, seed_seq,
                motif_set={k: v for k, v in NEURAL_MOTIFS.items() if k != "NRSF"},
                n_rounds=N_MOTIF_ROUNDS,
                verbose=True,
            )
            improvement = best_score - seed_score
            print(f"  Final score: {best_score:.4f}  (improvement: +{improvement:.4f})")
            all_motif_rows.append({
                "region":      region_name,
                "strategy":    "motif_insert",
                "seed_chrom":  row["chrom"],
                "seed_start":  row["start"],
                "seed_end":    row["end"],
                "seed_score":  round(seed_score, 4),
                "final_score": round(best_score, 4),
                "improvement": round(improvement, 4),
                "n_insertions": len(history),
                "insertions":  "; ".join(
                    f"{h['factor']}@{h['position']}" for h in history
                ),
                "gc_content":  round(gc_content(best_seq), 4),
                "sequence":    best_seq,
            })

    # ----------------------------------------------------------------
    # Save results
    # ----------------------------------------------------------------
    grad_path  = os.path.join(OUT_DIR, "gradient_designs.tsv")
    motif_path = os.path.join(OUT_DIR, "motif_designs.tsv")

    grad_df  = pd.DataFrame(all_gradient_rows).sort_values("cnn_score", ascending=False)
    motif_df = pd.DataFrame(all_motif_rows).sort_values("final_score", ascending=False)

    grad_df.to_csv(grad_path,  sep="\t", index=False)
    motif_df.to_csv(motif_path, sep="\t", index=False)

    print(f"\n{'='*60}")
    print("=== Design Summary ===")
    print(f"\nGradient designs saved to: {grad_path}")
    if len(grad_df):
        print(f"  Best score: {grad_df['cnn_score'].max():.4f}  "
              f"Mean: {grad_df['cnn_score'].mean():.4f}  "
              f"GC range: {grad_df['gc_content'].min():.2%}–{grad_df['gc_content'].max():.2%}")

    print(f"\nMotif-insertion designs saved to: {motif_path}")
    if len(motif_df):
        print(f"  Mean improvement: +{motif_df['improvement'].mean():.4f}")
        print(f"  Best final score: {motif_df['final_score'].max():.4f}")
        print(f"  Motifs inserted (example): {motif_df.iloc[0]['insertions']}")

    print(f"\nNext step: run python validate_designs.py (requires AlphaGenome API)")
    print("Or open data/processed/gradient_designs.tsv to inspect sequences.")


if __name__ == "__main__":
    main()

"""
Enhancer Identification Visualization

Produces a multi-panel figure per region:
  Panel 1  — Hi-C contact heatmap with TAD boundary lines
  Panel 2  — Insulation score track (50 kb window)
  Panel 3  — Enhancer candidate track (contact score per bin, coloured
              by enhancer_score when a trained model is available)

Usage
-----
    python visualize_enhancers.py

Requires
--------
- Mouse Micro-C .mcool at FILE_PATH
- Reference genome FASTA at FASTA_PATH (mm10; optional for contact-only mode)
- A trained EnhancerCNN at MODEL_PATH (optional; skips CNN scoring if missing)

If FASTA_PATH or MODEL_PATH is not found the script runs in contact-only mode,
which plots candidate bins ranked purely by contact enrichment.
"""

import os
import cooler
import cooltools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.gridspec import GridSpec

from src.tad_boundaries import (
    score_boundary_prominence,
    call_tad_intervals,
    multiscale_insulation,
    parse_coordinates,
    make_view_df,
)
from src.enhancer_id import extract_candidate_enhancers

# ---------------------------------------------------------------------------
# Configuration — update paths as needed
# ---------------------------------------------------------------------------

FILE_PATH  = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
FASTA_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mm10.fa"
MODEL_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/processed/enhancer_cnn.pt"

RESOLUTION        = 5000    # 5 kb
INSULATION_WINDOW = 50_000  # 50 kb

regions = {
    "Sox11_Chr12":   "chr12:26,000,000-28,000,000",
    "Mir9-2_Chr13":  "chr13:83,500,000-84,500,000",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _load_model():
    """Load EnhancerCNN weights if available; return None otherwise."""
    if not os.path.exists(MODEL_PATH):
        return None
    try:
        from src.enhancer_id import EnhancerCNN
        model = EnhancerCNN()
        model.load(MODEL_PATH)
        return model
    except Exception as e:
        print(f"  Note: could not load model from {MODEL_PATH}: {e}")
        return None


def _cnn_score_candidates(model, candidates: pd.DataFrame) -> pd.DataFrame:
    """Add enhancer_score column using the CNN if FASTA is available."""
    if model is None or not os.path.exists(FASTA_PATH):
        candidates['enhancer_score'] = np.nan
        return candidates
    try:
        from src.enhancer_id import fetch_sequences
        import torch
        seqs = fetch_sequences(candidates, FASTA_PATH, window_bp=500)
        model.model.eval()
        with torch.no_grad():
            X = torch.from_numpy(seqs).float()
            scores = model.model(X).squeeze(1).cpu().numpy()
        candidates['enhancer_score'] = scores.astype(np.float32)
    except Exception as e:
        print(f"  Note: CNN scoring skipped ({e})")
        candidates['enhancer_score'] = np.nan
    return candidates


# ---------------------------------------------------------------------------
# Main plotting function
# ---------------------------------------------------------------------------

def plot_enhancers(file_path: str, region_name: str, coordinates: str, model) -> None:
    print(f"Processing enhancers for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"  Error loading cooler: {e}")
        return

    chrom, start_bp, end_bp = parse_coordinates(coordinates)

    # --- 1. Contact matrix ---
    matrix = clr.matrix(balance=True).fetch(coordinates)
    n_bins = matrix.shape[0]

    # --- 2. Insulation + boundaries + TADs ---
    print("  Computing insulation and boundaries...")
    ins_table, windows = multiscale_insulation(
        clr, coordinates, window_sizes_bp=[INSULATION_WINDOW]
    )
    boundary_df = score_boundary_prominence(
        ins_table, INSULATION_WINDOW, prominence_thresholds=(0.2, 0.5)
    )
    tad_df = call_tad_intervals(boundary_df, clr, boundary_class_min='weak')

    ins_col = f'log2_insulation_score_{INSULATION_WINDOW}'
    ins_values = ins_table[ins_col].values if ins_col in ins_table.columns else np.zeros(n_bins)

    # --- 3. Candidate enhancers ---
    print("  Extracting candidate enhancers...")
    candidates = extract_candidate_enhancers(
        clr, coordinates, tad_df, boundary_df,
        contact_window_bp=50_000,
        boundary_exclusion_bp=10_000,
    )
    candidates = _cnn_score_candidates(model, candidates.copy())
    print(f"  Found {len(candidates)} candidate enhancer bins.")

    # --- 4. Build figure ---
    fig = plt.figure(figsize=(14, 11))
    gs = GridSpec(3, 1, height_ratios=[3, 1, 1], hspace=0.08, figure=fig)

    # Panel 1: Contact heatmap
    ax_heat = fig.add_subplot(gs[0])
    valid = matrix[np.isfinite(matrix)]
    vmin = np.percentile(valid, 2)  if len(valid) > 0 else 1e-4
    vmax = np.percentile(valid, 98) if len(valid) > 0 else 0.05
    vmin = max(vmin, 1e-5)
    im = ax_heat.imshow(
        matrix, cmap='RdYlBu_r',
        norm=mcolors.LogNorm(vmin=vmin, vmax=vmax),
        interpolation='none', aspect='auto',
    )
    ax_heat.set_title(f"Enhancer Candidates — {region_name}\n{coordinates}", fontsize=11)
    ax_heat.set_xticks([])
    ax_heat.set_ylabel("Genomic bins")
    plt.colorbar(im, ax=ax_heat, label="Contact frequency", pad=0.02, fraction=0.02)

    # Overlay strong boundary lines
    strong = boundary_df[boundary_df['boundary_class'] == 'strong']
    for _, row in strong.iterrows():
        pos = (row['start'] - start_bp) / RESOLUTION
        ax_heat.axvline(pos, color='white', lw=0.8, alpha=0.7)
        ax_heat.axhline(pos, color='white', lw=0.8, alpha=0.7)

    # Panel 2: Insulation score
    ax_ins = fig.add_subplot(gs[1])
    x_ins = np.arange(len(ins_values))
    ax_ins.plot(x_ins, ins_values, color='steelblue', lw=1.2, label=f'{INSULATION_WINDOW//1000}kb')
    ax_ins.fill_between(x_ins, ins_values, 0, where=(ins_values < 0),
                        color='steelblue', alpha=0.25, interpolate=True)
    ax_ins.axhline(0, color='grey', lw=0.5, ls='--')
    for _, row in strong.iterrows():
        pos = (row['start'] - start_bp) / RESOLUTION
        ax_ins.axvline(pos, color='red', lw=0.8, alpha=0.5)
    ax_ins.set_ylabel("Insulation\nscore", fontsize=9)
    ax_ins.set_xlim(0, n_bins)
    ax_ins.set_xticks([])
    ax_ins.legend(loc='lower right', fontsize='x-small')

    # Panel 3: Candidate enhancer track
    ax_enh = fig.add_subplot(gs[2])
    if len(candidates) > 0:
        bin_positions = ((candidates['start'] - start_bp) / RESOLUTION).values
        has_cnn = candidates['enhancer_score'].notna().any()

        if has_cnn:
            # Colour bars by CNN enhancer score
            norm = mcolors.Normalize(vmin=0.0, vmax=1.0)
            cmap = plt.cm.RdYlGn
            for i, (pos, score, cscore) in enumerate(zip(
                bin_positions,
                candidates['contact_score'].values,
                candidates['enhancer_score'].values,
            )):
                color = cmap(norm(cscore)) if np.isfinite(cscore) else 'grey'
                ax_enh.bar(pos, score, width=0.8, color=color, alpha=0.85, linewidth=0)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm.set_array([])
            plt.colorbar(sm, ax=ax_enh, label="CNN score", pad=0.02, fraction=0.02)
        else:
            # Contact-only mode
            ax_enh.bar(bin_positions, candidates['contact_score'].values,
                       width=0.8, color='darkorange', alpha=0.7, linewidth=0)

        # Mark bins called as enhancers (CNN score >= 0.5)
        if has_cnn:
            called = candidates[candidates['enhancer_score'] >= 0.5]
            if len(called) > 0:
                called_pos = ((called['start'] - start_bp) / RESOLUTION).values
                ax_enh.scatter(
                    called_pos,
                    called['contact_score'].values * 1.05,
                    marker='*', color='gold', edgecolors='black',
                    s=80, zorder=5, linewidth=0.5, label='enhancer (score≥0.5)',
                )
                ax_enh.legend(loc='upper right', fontsize='x-small')

    ax_enh.set_ylabel("Contact\nenrichment", fontsize=9)
    ax_enh.set_xlabel("Genomic bins")
    ax_enh.set_xlim(0, n_bins)

    # Save
    os.makedirs("media", exist_ok=True)
    output_path = f"media/{region_name}_enhancers.png"
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"  Saved → {output_path}")
    plt.close()


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if not os.path.exists(FILE_PATH):
        print(f"STOP: mcool file not found at {FILE_PATH}")
    else:
        if not os.path.exists(FASTA_PATH):
            print(f"Note: FASTA not found at {FASTA_PATH} — running in contact-only mode.")
        if not os.path.exists(MODEL_PATH):
            print(f"Note: no trained model at {MODEL_PATH} — running in contact-only mode.")

        model = _load_model()

        for name, coords in regions.items():
            plot_enhancers(FILE_PATH, name, coords, model)

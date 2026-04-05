import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os
import pandas as pd
from matplotlib.gridspec import GridSpec
from src.loaders import get_cooler

# --- CONFIGURATION ---
FILE_PATH = "mouse_microc.mcool"
RESOLUTION = 5000  # 5kb
REGION = "chr12:26,000,000-28,000,000" # Sox11 region
REGION_NAME = "Sox11_Locus"

# Deletion parameters (relative to region start)
# Let's delete a ~100kb region in the middle
DEL_START_BIN = 180 
DEL_END_BIN = 200   # 20 bins * 5kb = 100kb deletion
WINDOW_SIZE_BINS = 20 # 20*5kb = 100kb window for insulation

def compute_insulation_matrix(matrix, window_bins):
    """
    Computes insulation score for each bin in a matrix.
    Insulation(i) = average of contacts in a [W x W] diamond/square 
    centered on the diagonal at (i, i).
    """
    n = matrix.shape[0]
    insulation = np.full(n, np.nan)
    
    # Mirror matrix across diagonal just in case it's upper-triangular
    mat = np.nan_to_num(matrix)
    mat = (mat + mat.T) / 2.0
    
    for i in range(window_bins, n - window_bins):
        # The diamond is the region to the top-right of the diagonal point (i, i)
        # representing contacts stretching across bin i.
        # Specifically: matrix[i-W:i, i:i+W]
        square = mat[i - window_bins : i, i : i + window_bins]
        insulation[i] = np.nanmean(square)
    
    # Log2 normalization relative to the mean of the region
    avg_ins = np.nanmean(insulation)
    if avg_ins > 0:
        insulation = np.log2(insulation / avg_ins)
        
    return insulation

def apply_deletion(matrix, start_bin, end_bin):
    """
    Simulates a deletion by removing the rows and columns
    and shifting the rest of the matrix.
    """
    mat = matrix.copy()
    
    # Mask out the deleted region
    # A true deletion would actually remove the rows/cols.
    # We will do both: create a viz version (masked) and a calc version (sliced)
    mask = np.ones(mat.shape[0], dtype=bool)
    mask[start_bin:end_bin] = False
    
    # The 'edited' matrix is actual physical deletion (smaller size)
    edited_mat = mat[mask, :][:, mask]
    
    return edited_mat, mask

def run_prediction():
    print(f"Loading {FILE_PATH} at {RESOLUTION}bp...")
    try:
        clr = get_cooler(FILE_PATH, resolution=RESOLUTION)
    except Exception as e:
        print(f"Error: {e}")
        return

    print(f"Fetching region {REGION}...")
    matrix = clr.matrix(balance=True).fetch(REGION)
    n = matrix.shape[0]
    
    # 1. Baseline Insulation
    print(f"Computing baseline insulation (W={WINDOW_SIZE_BINS} bins)...")
    ins_ref = compute_insulation_matrix(matrix, WINDOW_SIZE_BINS)
    
    # 2. Simulate Deletion
    print(f"Simulating deletion of bins {DEL_START_BIN}-{DEL_END_BIN}...")
    edited_matrix, mask = apply_deletion(matrix, DEL_START_BIN, DEL_END_BIN)
    
    # 3. Post-Edit Insulation
    print("Computing edited insulation...")
    ins_edit_short = compute_insulation_matrix(edited_matrix, WINDOW_SIZE_BINS)
    
    # Maps edited bins back to original coordinate space for comparison
    ins_edit = np.full(n, np.nan)
    ins_edit[mask] = ins_edit_short
    
    # 4. Quantification
    delta = ins_edit - ins_ref
    
    # 5. Visualization
    print("Generating visualization...")
    fig = plt.figure(figsize=(15, 12))
    gs = GridSpec(4, 2, height_ratios=[4, 1, 1, 1], hspace=0.3)
    
    # Heatmap Ref
    ax_ref = fig.add_subplot(gs[0, 0])
    im1 = ax_ref.imshow(matrix, cmap='RdYlBu_r', norm=colors.LogNorm(vmin=0.001, vmax=0.05))
    ax_ref.set_title("Reference Matrix")
    
    # Heatmap Edit (Masked for visual alignment)
    ax_edit_v = fig.add_subplot(gs[0, 1])
    # For visualization, we'll show the same matrix but with a grey gap
    viz_mat = matrix.copy()
    viz_mat[DEL_START_BIN:DEL_END_BIN, :] = np.nan
    viz_mat[:, DEL_START_BIN:DEL_END_BIN] = np.nan
    im2 = ax_edit_v.imshow(viz_mat, cmap='RdYlBu_r', norm=colors.LogNorm(vmin=0.001, vmax=0.05))
    ax_edit_v.set_title(f"Simulated Deletion (Bins {DEL_START_BIN}-{DEL_END_BIN})")
    
    # Insulation Tracks
    ax_ins = fig.add_subplot(gs[1, :])
    x = np.arange(n)
    ax_ins.plot(x, ins_ref, label="Reference Insulation", color='black', alpha=0.7)
    ax_ins.plot(x, ins_edit, label="Post-Deletion Insulation", color='red', lw=2)
    ax_ins.axvspan(DEL_START_BIN, DEL_END_BIN, color='grey', alpha=0.3, label="Deleted Region")
    ax_ins.set_ylabel("log2 Insulation")
    ax_ins.legend(loc='upper right', fontsize='small')
    ax_ins.set_title("Insulation Score Prediction")
    
    # Delta Track
    ax_delta = fig.add_subplot(gs[2, :])
    ax_delta.fill_between(x, 0, delta, where=(delta > 0), color='blue', alpha=0.5, label="Insulation Increase")
    ax_delta.fill_between(x, 0, delta, where=(delta < 0), color='orange', alpha=0.5, label="Insulation Decrease")
    ax_delta.axhline(0, color='black', lw=0.5)
    ax_delta.set_ylabel("Delta log2")
    ax_delta.legend(loc='upper right', fontsize='small')
    ax_delta.set_title("Structural Impact (Delta Insulation)")
    
    # Distribution of Delta
    ax_dist = fig.add_subplot(gs[3, :])
    valid_delta = delta[np.isfinite(delta)]
    ax_dist.hist(valid_delta, bins=50, color='purple', alpha=0.6)
    ax_dist.set_title("Distribution of Delta Insulation Across Locus")
    ax_dist.set_xlabel("Change in log2 Insulation")
    
    # Save
    os.makedirs("media", exist_ok=True)
    out_path = "media/insulation_delta_prediction.png"
    plt.savefig(out_path, dpi=200, bbox_inches='tight')
    plt.close()
    print(f"Success! Prediction saved to {out_path}")

if __name__ == "__main__":
    run_prediction()

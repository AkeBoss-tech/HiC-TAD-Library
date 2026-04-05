import os
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap

def plot_hic_triangle(ref_map: np.ndarray, mut_map: np.ndarray, output_path: str, title: str = "3D Chromatin Contact Map (Ref vs Mut)"):
    """
    Generate a triangular Hi-C heatmap for Ref vs Mut and their delta.

    Args:
        ref_map: 2D numpy array for reference contact frequencies.
        mut_map: 2D numpy array for mutant contact frequencies.
        output_path: Path to save the resulting image.
        title: Title of the plot.
    """
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))

    # We rotate the matrix to make it triangular for standard Hi-C view.
    # To keep it simple, we just plot the 2D matrix directly but use upper triangle.
    ref_tri = np.triu(ref_map)
    mut_tri = np.triu(mut_map)
    delta_tri = mut_tri - ref_tri

    cmap_hic = "Reds"
    cmap_delta = "coolwarm"

    im0 = axes[0].imshow(ref_tri, cmap=cmap_hic, interpolation='nearest')
    axes[0].set_title("Reference")
    plt.colorbar(im0, ax=axes[0], fraction=0.046, pad=0.04)

    im1 = axes[1].imshow(mut_tri, cmap=cmap_hic, interpolation='nearest')
    axes[1].set_title("Mutant")
    plt.colorbar(im1, ax=axes[1], fraction=0.046, pad=0.04)

    vmax = np.max(np.abs(delta_tri))
    im2 = axes[2].imshow(delta_tri, cmap=cmap_delta, vmin=-vmax, vmax=vmax, interpolation='nearest')
    axes[2].set_title("Delta (Mut - Ref)")
    plt.colorbar(im2, ax=axes[2], fraction=0.046, pad=0.04)

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path

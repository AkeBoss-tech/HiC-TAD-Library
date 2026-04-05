import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate


def _to_triangle(matrix: np.ndarray) -> np.ndarray:
    """
    Rotate an upper-triangle contact matrix 45° to produce the standard
    triangular Hi-C view (diagonal becomes the base).
    """
    tri = np.triu(matrix).astype(float)
    tri[tri == 0] = np.nan
    rotated = rotate(tri, angle=45, reshape=True, cval=np.nan, order=1)
    h = rotated.shape[0]
    return rotated[h // 2:, :]


def plot_hic_triangle(ref_map: np.ndarray, mut_map: np.ndarray, output_path: str,
                      title: str = "3D Chromatin Contact Map (Ref vs Mut)"):
    """
    Generate triangular Hi-C heatmaps for Ref, Mut, and their delta.

    The matrix is rotated 45° so the diagonal sits at the bottom edge,
    matching the standard genomics convention.

    Args:
        ref_map: 2D numpy array for reference contact frequencies.
        mut_map: 2D numpy array for mutant contact frequencies.
        output_path: Path to save the resulting image.
        title: Title of the plot.
    """
    ref_tri = _to_triangle(ref_map)
    mut_tri = _to_triangle(mut_map)

    delta_raw = np.triu(mut_map - ref_map).astype(float)
    delta_raw[np.triu(ref_map) == 0] = np.nan
    delta_tri = rotate(delta_raw, angle=45, reshape=True, cval=np.nan, order=1)
    h = delta_tri.shape[0]
    delta_tri = delta_tri[h // 2:, :]

    fig, axes = plt.subplots(3, 1, figsize=(10, 9))

    cmap_hic = "Reds"
    cmap_delta = "coolwarm"

    for ax, mat, label in zip(axes[:2], [ref_tri, mut_tri], ["Reference", "Mutant"]):
        im = ax.imshow(mat, cmap=cmap_hic, interpolation="nearest", aspect="auto",
                       origin="upper")
        ax.set_title(label)
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)

    vmax = np.nanmax(np.abs(delta_tri)) if not np.all(np.isnan(delta_tri)) else 1.0
    im2 = axes[2].imshow(delta_tri, cmap=cmap_delta, vmin=-vmax, vmax=vmax,
                         interpolation="nearest", aspect="auto", origin="upper")
    axes[2].set_title("Delta (Mut − Ref)")
    axes[2].set_yticks([])
    plt.colorbar(im2, ax=axes[2], fraction=0.02, pad=0.02)

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path

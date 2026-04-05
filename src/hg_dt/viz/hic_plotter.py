import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import rotate


def _to_triangle(matrix: np.ndarray) -> np.ndarray:
    """
    Rotate an upper-triangle contact matrix 45° CCW to produce the standard
    triangular Hi-C view (diagonal at the bottom, long-range contacts at top).

    After a 45° CCW rotation, the upper-triangle data lands in the TOP half of
    the rotated image. Taking [:h//2] gives us the triangle with the diagonal
    at the bottom row and long-range contacts ascending toward the top.
    """
    tri = np.triu(matrix).astype(float)
    # Only mask the strict lower triangle (below diagonal) — zeros are valid data
    rows, cols = np.tril_indices(tri.shape[0], k=-1)
    tri[rows, cols] = np.nan
    rotated = rotate(tri, angle=45, reshape=True, cval=np.nan, order=1)
    h = rotated.shape[0]
    # TOP half contains the upper-triangle data after 45° CCW rotation
    return rotated[:h // 2, :]


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
    # Mask strict lower triangle only (not zeros)
    rows, cols = np.tril_indices(delta_raw.shape[0], k=-1)
    delta_raw[rows, cols] = np.nan
    delta_tri = rotate(delta_raw, angle=45, reshape=True, cval=np.nan, order=1)
    h = delta_tri.shape[0]
    delta_tri = delta_tri[:h // 2, :]   # top half contains upper-triangle data

    fig, axes = plt.subplots(3, 1, figsize=(10, 9))

    cmap_hic   = "YlOrRd"   # yellow-orange-red, perceptually uniform for contacts
    cmap_delta = "RdBu_r"   # red=gained contacts, blue=lost (intuitive for delta)

    for ax, mat, label in zip(axes[:2], [ref_tri, mut_tri], ["Reference", "Mutant"]):
        vmax_hic = np.nanmax(mat) if not np.all(np.isnan(mat)) else 1.0
        im = ax.imshow(mat, cmap=cmap_hic, vmin=0, vmax=vmax_hic,
                       interpolation="nearest", aspect="auto", origin="upper")
        ax.set_title(label)
        ax.set_yticks([])
        plt.colorbar(im, ax=ax, fraction=0.02, pad=0.02)

    vmax = np.nanmax(np.abs(delta_tri)) if not np.all(np.isnan(delta_tri)) else 1.0
    im2 = axes[2].imshow(delta_tri, cmap=cmap_delta, vmin=-vmax, vmax=vmax,
                         interpolation="nearest", aspect="auto", origin="upper")
    axes[2].set_title("Delta (Mut − Ref)  ·  Red = gained contacts  ·  Blue = lost")
    axes[2].set_yticks([])
    plt.colorbar(im2, ax=axes[2], fraction=0.02, pad=0.02)

    plt.suptitle(title)
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path

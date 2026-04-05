import os
import matplotlib.pyplot as plt
import numpy as np

def plot_tracks(ref_track: np.ndarray, mut_track: np.ndarray, output_path: str, title: str = "1D Genome Tracks (Ref vs. Mut)"):
    """
    Generate a 1D browser plot comparing reference and mutant tracks.
    Uses matplotlib as a lightweight alternative to a full pyGenomeTracks setup
    for generating the 1D delta view (accessibility / expression).

    Args:
        ref_track: 1D numpy array of reference signal.
        mut_track: 1D numpy array of mutant signal.
        output_path: Path to save the resulting image.
        title: Title of the plot.
    """
    fig, axes = plt.subplots(3, 1, figsize=(10, 6), sharex=True)

    # Ref track
    axes[0].fill_between(range(len(ref_track)), ref_track, color="blue", alpha=0.7)
    axes[0].set_ylabel("Reference")
    axes[0].set_title(title)

    # Mut track
    axes[1].fill_between(range(len(mut_track)), mut_track, color="orange", alpha=0.7)
    axes[1].set_ylabel("Mutant")

    # Delta track
    delta = mut_track - ref_track
    axes[2].fill_between(range(len(delta)), delta, where=(delta >= 0), color="red", alpha=0.7, label="Gain")
    axes[2].fill_between(range(len(delta)), delta, where=(delta < 0), color="green", alpha=0.7, label="Loss")
    axes[2].set_ylabel("Delta (Mut - Ref)")
    axes[2].set_xlabel("Genomic Position (bins)")
    axes[2].legend()

    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path

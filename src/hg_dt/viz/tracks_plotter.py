import os
import matplotlib.pyplot as plt
import numpy as np

def plot_tracks(ref_track: np.ndarray, mut_track: np.ndarray, output_path: str, title: str = "1D Genome Tracks (Ref vs. Mut)"):
    """
    Generate a 1D browser plot comparing reference and mutant tracks.
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

def plot_multi_tracks(tracks_data: dict, output_path: str, title: str = "AlphaGenome Multi-Scale Predictions"):
    """
    Generate a vertically-stacked genome browser plot for multiple signals.
    
    Args:
        tracks_data: Dictionary mapping track name (e.g., 'ATAC', 'CTCF') to 
                    a tuple (ref_array, mut_array).
        output_path: Path to save the resulting image.
        title: Title of the entire figure.
    """
    n_tracks = len(tracks_data)
    fig, axes = plt.subplots(n_tracks, 1, figsize=(10, 3 * n_tracks), sharex=True, squeeze=False)
    
    for i, (track_name, (ref, mut)) in enumerate(tracks_data.items()):
        ax = axes[i, 0]
        # Reference in Grey/Base
        ax.fill_between(range(len(ref)), ref, color="grey", alpha=0.3, label="Reference")
        
        # Mutant in color matching biological convention
        color_map = {
            'ATAC': 'darkred',
            'CTCF': 'blue',
            'H3K27ac': 'orange',
            'H3K4me3': 'green',
            'CAGE': 'purple'
        }
        color = color_map.get(track_name, "blue")
        ax.plot(range(len(mut)), mut, color=color, linewidth=1.5, label="Mutant Prediction")
        
        # Delta indicators (gain/loss)
        delta = mut - ref
        ax.fill_between(range(len(delta)), mut, ref, where=(delta < 0), color="green", alpha=0.5, label="Loss")
        ax.fill_between(range(len(delta)), mut, ref, where=(delta > 0), color="red", alpha=0.5, label="Gain")
        
        ax.set_ylabel(f"{track_name}\n(Signal)")
        ax.legend(loc='upper right', fontsize='x-small')
        if i == 0:
            ax.set_title(title)
            
    axes[-1, 0].set_xlabel("Genomic Position (bins)")
    plt.tight_layout()
    plt.savefig(output_path, dpi=150)
    plt.close()
    return output_path

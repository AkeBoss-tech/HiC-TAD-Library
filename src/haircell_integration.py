"""
Integration functions for combining observed Hi-C data with AlphaGenome
predicted regulatory signals for hair cell differentiation analysis.

Provides:
- Insulation scoring on raw contact matrices (observed or predicted)
- Boundary detection from insulation arrays
- Concordance analysis between observed and predicted boundaries
- Virtual 4C extraction
- Visualization helpers for differential contacts, boundary validation,
  virtual 4C comparison, and ISM variant scores
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.signal import find_peaks
from typing import Tuple, Optional


# ---------------------------------------------------------------------------
# Analysis functions
# ---------------------------------------------------------------------------

def compute_insulation_from_matrix(
    contact_matrix: np.ndarray,
    window_bins: int = 10,
) -> np.ndarray:
    """
    Compute insulation score from a square contact matrix.

    Uses the diamond-window approach: for each diagonal position *i*, average
    the contacts in the window_bins x window_bins sub-matrix that straddles the
    diagonal.  Low values indicate insulation (boundary).

    Works on both observed (cooler) and predicted (AlphaGenome) matrices.

    Parameters
    ----------
    contact_matrix : np.ndarray
        Square contact frequency matrix.
    window_bins : int
        Half-size of the diamond window in bins.

    Returns
    -------
    np.ndarray
        Log2 insulation scores (NaN at edges where window doesn't fit).
    """
    n = contact_matrix.shape[0]
    insulation = np.full(n, np.nan)

    for i in range(window_bins, n - window_bins):
        diamond = contact_matrix[i - window_bins:i, i:i + window_bins]
        score = np.nanmean(diamond)
        if score > 0:
            insulation[i] = score

    mean_ins = np.nanmean(insulation)
    if mean_ins > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            insulation = np.log2(insulation / mean_ins)

    return insulation


def detect_boundaries_from_insulation(
    insulation: np.ndarray,
    prominence_thresholds: Tuple[float, float] = (0.2, 0.5),
    min_distance_bins: int = 3,
) -> pd.DataFrame:
    """
    Detect TAD boundaries from a 1-D insulation score array.

    Parameters
    ----------
    insulation : np.ndarray
        Log2 insulation scores (output of compute_insulation_from_matrix).
    prominence_thresholds : (float, float)
        (weak_min, strong_min) prominence cutoffs.
    min_distance_bins : int
        Minimum separation between detected minima.

    Returns
    -------
    pd.DataFrame
        Columns: bin_index, prominence, boundary_class.
    """
    ins = insulation.copy()
    ins[np.isnan(ins)] = 0.0

    peaks, props = find_peaks(-ins, distance=min_distance_bins, prominence=0)

    weak_min, strong_min = prominence_thresholds
    results = []
    for peak, prom in zip(peaks, props['prominences']):
        if prom >= strong_min:
            cls = 'strong'
        elif prom >= weak_min:
            cls = 'weak'
        else:
            cls = 'sub_threshold'
        results.append({
            'bin_index': peak,
            'prominence': prom,
            'boundary_class': cls,
        })

    if results:
        return pd.DataFrame(results)
    return pd.DataFrame(columns=['bin_index', 'prominence', 'boundary_class'])


def differential_contact_map(
    map_a: np.ndarray,
    map_b: np.ndarray,
) -> np.ndarray:
    """Element-wise difference of two contact maps (A - B), handling NaNs."""
    with np.errstate(invalid='ignore'):
        return map_a - map_b


def boundary_concordance(
    observed_bins: pd.DataFrame,
    predicted_bins: pd.DataFrame,
    tolerance_bins: int = 3,
) -> pd.DataFrame:
    """
    Assess concordance between observed and predicted boundary positions.

    A boundary is *concordant* (sequence-encoded) if a predicted boundary
    falls within *tolerance_bins* of the observed position.

    Parameters
    ----------
    observed_bins : pd.DataFrame
        Must contain 'bin_index', 'boundary_class', 'prominence'.
    predicted_bins : pd.DataFrame
        Must contain 'bin_index', 'boundary_class', 'prominence'.
    tolerance_bins : int
        Maximum distance (in bins) to count as a match.

    Returns
    -------
    pd.DataFrame
        One row per observed boundary with match status and interpretation.
    """
    results = []

    for _, obs in observed_bins.iterrows():
        obs_pos = int(obs['bin_index'])

        if len(predicted_bins) > 0:
            distances = np.abs(predicted_bins['bin_index'].values - obs_pos)
            nearest_idx = np.argmin(distances)
            nearest_dist = int(distances[nearest_idx])
            nearest_pred = predicted_bins.iloc[nearest_idx]
            matched = nearest_dist <= tolerance_bins
        else:
            nearest_dist = -1
            nearest_pred = None
            matched = False

        results.append({
            'observed_bin': obs_pos,
            'observed_class': obs['boundary_class'],
            'observed_prominence': obs['prominence'],
            'nearest_predicted_bin': int(nearest_pred['bin_index']) if matched else None,
            'predicted_class': nearest_pred['boundary_class'] if matched else None,
            'distance_bins': nearest_dist,
            'is_concordant': matched,
            'interpretation': (
                'sequence-encoded' if matched else 'epigenetically regulated'
            ),
        })

    return pd.DataFrame(results)


def extract_virtual_4c(
    contact_matrix: np.ndarray,
    viewpoint_bin: int,
) -> np.ndarray:
    """Extract one row of a contact matrix as a virtual-4C profile."""
    return contact_matrix[viewpoint_bin, :]


# ---------------------------------------------------------------------------
# Visualization
# ---------------------------------------------------------------------------

def plot_differential_contacts(
    diff_map: np.ndarray,
    interval_str: str,
    title: str = "Differential Contact Map (IHC - OHC)",
    output_path: Optional[str] = None,
):
    """Heatmap of a differential contact map with diverging colormap."""
    fig, ax = plt.subplots(figsize=(8, 8))

    vmax = np.nanpercentile(np.abs(diff_map), 95)
    if vmax == 0:
        vmax = 1.0

    im = ax.imshow(
        diff_map, cmap='RdBu_r', interpolation='none',
        vmin=-vmax, vmax=vmax,
    )
    ax.set_title(f"{title}\n{interval_str}")
    plt.colorbar(im, ax=ax, label="Contact Difference (IHC - OHC)")

    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_path}")
    plt.close()


def plot_boundary_validation(
    predicted_insulation: np.ndarray,
    concordance: pd.DataFrame,
    interval_str: str,
    resolution: int,
    title: str = "Boundary Validation",
    output_path: Optional[str] = None,
):
    """Predicted insulation track annotated with concordance markers."""
    fig, axes = plt.subplots(
        2, 1, figsize=(14, 6),
        gridspec_kw={'height_ratios': [3, 1]}, sharex=True,
    )

    x = np.arange(len(predicted_insulation))

    # Insulation track
    axes[0].plot(x, predicted_insulation, color='black', lw=1)
    axes[0].fill_between(x, predicted_insulation, 0, alpha=0.3, color='steelblue')
    axes[0].set_ylabel("Predicted Insulation\n(AlphaGenome)")
    axes[0].set_title(f"{title}\n{interval_str}")

    # Mark boundaries
    for _, row in concordance.iterrows():
        color = 'green' if row['is_concordant'] else 'red'
        axes[0].axvline(row['observed_bin'], color=color, lw=1.5, alpha=0.6, ls='--')
        axes[1].axvline(row['observed_bin'], color=color, lw=2, alpha=0.7)

    axes[1].set_ylabel("Boundary\nType")
    axes[1].set_xlabel(f"Position (bins, {resolution}bp)")

    legend_elements = [
        Line2D([0], [0], color='green', lw=2, label='Sequence-encoded (conserved)'),
        Line2D([0], [0], color='red', lw=2, label='Epigenetically regulated'),
    ]
    axes[1].legend(handles=legend_elements, loc='upper right', fontsize='small')

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_path}")
    plt.close()


def plot_virtual_4c_comparison(
    v4c_a: np.ndarray,
    v4c_b: np.ndarray,
    viewpoint_bin: int,
    interval_str: str,
    label_a: str = "IHC",
    label_b: str = "OHC",
    gene_name: str = "",
    output_path: Optional[str] = None,
):
    """Overlay and difference of two virtual-4C profiles."""
    fig, axes = plt.subplots(2, 1, figsize=(14, 6), sharex=True)

    min_len = min(len(v4c_a), len(v4c_b))
    x = np.arange(min_len)

    # Overlay
    axes[0].fill_between(x, v4c_a[:min_len], alpha=0.5, color='steelblue', label=label_a)
    axes[0].fill_between(x, v4c_b[:min_len], alpha=0.5, color='coral', label=label_b)
    axes[0].axvline(viewpoint_bin, color='black', ls='--', lw=1, label='Viewpoint')
    axes[0].set_ylabel("Contact Frequency")
    axes[0].legend(loc='upper right')
    axes[0].set_title(f"Virtual 4C: {gene_name} promoter\n{interval_str}")

    # Difference
    diff = v4c_a[:min_len] - v4c_b[:min_len]
    axes[1].fill_between(
        x, diff, 0, where=(diff > 0),
        color='steelblue', alpha=0.5, label=f'{label_a} > {label_b}',
    )
    axes[1].fill_between(
        x, diff, 0, where=(diff < 0),
        color='coral', alpha=0.5, label=f'{label_b} > {label_a}',
    )
    axes[1].axvline(viewpoint_bin, color='black', ls='--', lw=1)
    axes[1].axhline(0, color='grey', lw=0.5)
    axes[1].set_ylabel(f"Difference\n({label_a} - {label_b})")
    axes[1].set_xlabel("Position (bins)")
    axes[1].legend(loc='upper right')

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_path}")
    plt.close()


def plot_ism_scores(
    scores_df: pd.DataFrame,
    locus_name: str,
    output_path: Optional[str] = None,
):
    """Bar chart of ISM variant scores ranked by TAD dissolution impact."""
    if len(scores_df) == 0:
        return

    sorted_df = scores_df.sort_values(
        'tad_dissolution_score', ascending=False,
    ).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(12, 5))
    bars = ax.bar(
        range(len(sorted_df)),
        sorted_df['tad_dissolution_score'].values,
        color='steelblue', alpha=0.7,
    )

    # Highlight the top 5 most disruptive variants
    for i in range(min(5, len(sorted_df))):
        bars[i].set_color('firebrick')
        bars[i].set_alpha(1.0)
        row = sorted_df.iloc[i]
        ax.annotate(
            f"{row['ref']}>{row['alt']}",
            (i, row['tad_dissolution_score']),
            fontsize=7, ha='center', va='bottom', rotation=45,
        )

    ax.set_xlabel("Variant (ranked by effect size)")
    ax.set_ylabel("TAD Dissolution Score\n(mean |insulation change|)")
    ax.set_title(f"In Silico Mutagenesis: {locus_name} Boundary CTCF Sites")

    plt.tight_layout()
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"  Saved to {output_path}")
    plt.close()

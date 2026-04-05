import numpy as np
from typing import Dict, Tuple, List

def accessibility_delta(ref_atac: np.ndarray, mut_atac: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """
    Compute log2(Mut/Ref) for accessibility across the window.
    (Note: Log2(Mut/Ref) so negative means loss of accessibility in mutant)

    Args:
        ref_atac: The predicted reference ATAC-seq track.
        mut_atac: The predicted mutant ATAC-seq track.
        pseudocount: Small value to avoid division by zero or log of zero.

    Returns:
        np.ndarray containing the log2 fold change.
    """
    ref = ref_atac + pseudocount
    mut = mut_atac + pseudocount
    return np.log2(mut / ref)

def contact_delta(ref_map: np.ndarray, mut_map: np.ndarray) -> np.ndarray:
    """
    Compute difference heatmap for predicted loops.

    Args:
        ref_map: Predicted reference contact map (2D numpy array).
        mut_map: Predicted mutant contact map (2D numpy array).

    Returns:
        np.ndarray containing the difference (mut - ref).
    """
    return mut_map - ref_map

def expression_delta(ref_expr: np.ndarray, mut_expr: np.ndarray, pseudocount: float = 1e-6) -> np.ndarray:
    """
    Compute predicted RNA abundance change for target promoters.

    Args:
        ref_expr: Predicted reference RNA-seq or CAGE track.
        mut_expr: Predicted mutant RNA-seq or CAGE track.
        pseudocount: Small value to avoid division by zero or log of zero.

    Returns:
        np.ndarray containing the log2 fold change.
    """
    ref = ref_expr + pseudocount
    mut = mut_expr + pseudocount
    return np.log2(mut / ref)

def find_silenced_elements(ref_atac: np.ndarray, mut_atac: np.ndarray, threshold: float = 0.5) -> List[int]:
    """
    Find elements where predicted accessibility drops by >50%.

    Args:
        ref_atac: Predicted reference ATAC-seq track.
        mut_atac: Predicted mutant ATAC-seq track.
        threshold: The fractional drop threshold (e.g., 0.5 for >50% drop).

    Returns:
        List of indices where the accessibility dropped by more than the threshold.
    """
    # >50% drop means mut < ref * (1 - threshold)
    # Ensure we only consider regions where ref accessibility is meaningful (> 0)
    significant_ref = ref_atac > 1e-3
    silenced = (mut_atac < ref_atac * (1.0 - threshold)) & significant_ref
    return np.where(silenced)[0].tolist()

def find_distal_loops(contact_map: np.ndarray, min_dist: int = 5, min_freq: float = 0.1) -> List[Tuple[int, int]]:
    """
    Identify promoter-enhancer pairs (loops) in the contact map.

    Args:
        contact_map: 2D numpy array of contact frequencies.
        min_dist: Minimum distance (in bins) from the diagonal to be considered a distal loop.
        min_freq: Minimum contact frequency threshold.

    Returns:
        List of tuples (i, j) representing the bins involved in the loop.
    """
    loops = []
    rows, cols = contact_map.shape
    for i in range(rows):
        for j in range(i + min_dist, cols):
            if contact_map[i, j] > min_freq:
                loops.append((i, j))
    return loops

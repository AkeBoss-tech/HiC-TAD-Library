"""
Lightweight 3-D polymer simulation from Hi-C contact matrices.

Converts a balanced contact matrix into harmonic spring restraints and
runs overdamped Langevin dynamics to produce bead coordinates.  No
external simulation engine (polychrom, OpenMM) required — only numpy
and scipy.

Typical workflow
----------------
>>> from src.polymer_sim import contact_matrix_to_restraints, simulate_polymer
>>> matrix = clr.matrix(balance=True).fetch("chr12:26,000,000-28,000,000")
>>> restraints = contact_matrix_to_restraints(matrix)
>>> coords = simulate_polymer(len(matrix), restraints)
"""

import numpy as np
from typing import List, Tuple, Optional


# ---------------------------------------------------------------------------
# 1. Convert contact matrix → spring restraints
# ---------------------------------------------------------------------------

def contact_matrix_to_restraints(
    matrix: np.ndarray,
    contact_threshold_quantile: float = 0.70,
    k_min: float = 0.1,
    k_max: float = 5.0,
) -> List[Tuple[int, int, float, float]]:
    """
    Derive harmonic spring restraints from a balanced Hi-C contact matrix.

    Each entry above the threshold becomes a restraint (i, j, rest_length, k)
    where *k* scales with log-contact-frequency and *rest_length* decreases
    for stronger contacts.

    Parameters
    ----------
    matrix : np.ndarray
        Square, balanced contact matrix (NaN-safe).
    contact_threshold_quantile : float
        Only contacts above this quantile of the off-diagonal values
        are converted to restraints.  Lower = more springs, slower sim.
    k_min, k_max : float
        Minimum and maximum spring constants.

    Returns
    -------
    list of (i, j, rest_length, k)
        Each tuple describes one harmonic restraint between beads i and j.
    """
    mat = np.nan_to_num(matrix, nan=0.0).copy()
    np.fill_diagonal(mat, 0.0)
    # Symmetrise
    mat = (mat + mat.T) / 2.0

    # Only upper triangle (avoid double-counting)
    triu_i, triu_j = np.triu_indices_from(mat, k=2)  # skip ±1 diagonal (backbone)
    values = mat[triu_i, triu_j]
    positive = values > 0
    triu_i, triu_j, values = triu_i[positive], triu_j[positive], values[positive]

    if len(values) == 0:
        return []

    threshold = np.quantile(values, contact_threshold_quantile)
    keep = values >= threshold
    triu_i, triu_j, values = triu_i[keep], triu_j[keep], values[keep]

    if len(values) == 0:
        return []

    log_vals = np.log1p(values)
    lo, hi = log_vals.min(), log_vals.max()
    if hi == lo:
        norm = np.ones_like(log_vals) * 0.5
    else:
        norm = (log_vals - lo) / (hi - lo)

    restraints: List[Tuple[int, int, float, float]] = []
    for idx in range(len(triu_i)):
        i, j = int(triu_i[idx]), int(triu_j[idx])
        t = norm[idx]
        k = k_min + t * (k_max - k_min)
        # Stronger contact → shorter rest length (range 1.0 – 3.0 bead-units)
        rest_length = 3.0 - 2.0 * t
        restraints.append((i, j, rest_length, k))

    return restraints


def insulation_to_backbone_stiffness(
    insulation_scores: np.ndarray,
    k_soft: float = 1.0,
    k_stiff: float = 10.0,
) -> np.ndarray:
    """
    Map per-bin insulation scores to backbone spring constants.

    Lower insulation (deeper valley = stronger boundary) → stiffer spring
    between consecutive beads, mimicking the observation that TAD boundaries
    constrain spatial proximity.

    Parameters
    ----------
    insulation_scores : np.ndarray
        log2 insulation scores for N bins. NaNs are replaced with 0.
    k_soft, k_stiff : float
        Range of backbone spring constants.

    Returns
    -------
    np.ndarray
        Array of length N-1 with spring constant for each consecutive pair.
    """
    scores = np.nan_to_num(insulation_scores, nan=0.0)
    # Normalise: most-negative → stiffest, most-positive → softest
    lo, hi = scores.min(), scores.max()
    if hi == lo:
        return np.full(len(scores) - 1, (k_soft + k_stiff) / 2)
    norm = (scores - lo) / (hi - lo)  # 0 = low insulation, 1 = high
    # Pair-wise: average of neighbors
    pair_norm = (norm[:-1] + norm[1:]) / 2.0
    # Invert: low insulation → high stiffness
    k = k_stiff - pair_norm * (k_stiff - k_soft)
    return k


# ---------------------------------------------------------------------------
# 2. Simulation engine
# ---------------------------------------------------------------------------

def _init_coords(n_beads: int, rng: np.random.Generator) -> np.ndarray:
    """Initialise beads along a random-walk backbone."""
    coords = np.zeros((n_beads, 3))
    for i in range(1, n_beads):
        step = rng.standard_normal(3)
        step /= np.linalg.norm(step) + 1e-12
        coords[i] = coords[i - 1] + step
    return coords


def simulate_polymer(
    n_beads: int,
    restraints: List[Tuple[int, int, float, float]],
    backbone_k: Optional[np.ndarray] = None,
    backbone_rest: float = 1.0,
    n_steps: int = 5000,
    dt: float = 0.005,
    friction: float = 1.0,
    temperature: float = 1.0,
    seed: int = 42,
) -> np.ndarray:
    """
    Run overdamped Langevin dynamics for a bead-spring polymer.

    Forces
    ------
    * Backbone: harmonic springs between consecutive beads.
    * Hi-C restraints: harmonic springs from ``contact_matrix_to_restraints``.
    * Excluded-volume: soft repulsive potential (truncated LJ-like).
    * Thermal noise: Gaussian kicks scaled by temperature.

    Parameters
    ----------
    n_beads : int
        Number of beads (= number of genomic bins in the region).
    restraints : list of (i, j, rest_length, k)
        Non-backbone spring restraints derived from Hi-C contacts.
    backbone_k : np.ndarray or None
        Per-pair backbone stiffness (length n_beads-1).  If None, uses a
        uniform stiffness of 5.0 for every pair.
    backbone_rest : float
        Rest length of backbone springs.
    n_steps : int
        Number of integration steps.
    dt : float
        Time step.
    friction : float
        Damping coefficient.
    temperature : float
        Thermal energy scale (kT).
    seed : int
        Random seed for reproducibility.

    Returns
    -------
    np.ndarray
        Final bead positions, shape (n_beads, 3).
    """
    rng = np.random.default_rng(seed)
    coords = _init_coords(n_beads, rng)

    if backbone_k is None:
        backbone_k = np.full(n_beads - 1, 5.0)

    noise_scale = np.sqrt(2.0 * temperature * dt / friction)

    # Pre-convert restraints to arrays for speed
    if restraints:
        r_i = np.array([r[0] for r in restraints], dtype=int)
        r_j = np.array([r[1] for r in restraints], dtype=int)
        r_d0 = np.array([r[2] for r in restraints], dtype=float)
        r_k = np.array([r[3] for r in restraints], dtype=float)
    else:
        r_i = np.array([], dtype=int)
        r_j = np.array([], dtype=int)
        r_d0 = np.array([], dtype=float)
        r_k = np.array([], dtype=float)

    for _ in range(n_steps):
        forces = np.zeros_like(coords)

        # --- Backbone springs ---
        diff = coords[1:] - coords[:-1]                    # (n-1, 3)
        dists = np.linalg.norm(diff, axis=1, keepdims=True) + 1e-12
        unit = diff / dists
        stretch = dists.squeeze() - backbone_rest
        f_mag = backbone_k * stretch
        f_backbone = (f_mag[:, None] * unit)
        forces[:-1] += f_backbone
        forces[1:] -= f_backbone

        # --- Hi-C restraints ---
        if len(r_i) > 0:
            diff_r = coords[r_j] - coords[r_i]
            dist_r = np.linalg.norm(diff_r, axis=1, keepdims=True) + 1e-12
            unit_r = diff_r / dist_r
            stretch_r = dist_r.squeeze() - r_d0
            f_r = (r_k * stretch_r)[:, None] * unit_r
            np.add.at(forces, r_i, f_r)
            np.add.at(forces, r_j, -f_r)

        # --- Soft excluded-volume (only nearby beads) ---
        # Brute-force is fine for N < ~2000 beads
        if n_beads <= 2000:
            for i in range(n_beads):
                d = coords[i + 2:] - coords[i]
                dn = np.linalg.norm(d, axis=1, keepdims=True) + 1e-12
                overlap = dn.squeeze() < 0.8
                if np.any(overlap):
                    rep = d[overlap] / dn[overlap]
                    # Explicitly reshape to ensure proper broadcasting: (m, 1) shape
                    dist_diff = (0.8 - dn[overlap]).reshape(-1, 1)
                    rep_force = rep * dist_diff * 10.0
                    forces[i] -= rep_force.sum(axis=0).ravel()
                    forces[i + 2:][overlap] += rep_force

        # --- Integration (overdamped: v = F / gamma + noise) ---
        coords += (forces / friction) * dt + rng.standard_normal(coords.shape) * noise_scale

    return coords


# ---------------------------------------------------------------------------
# 3. Convenience wrapper: cooler region → 3-D coordinates
# ---------------------------------------------------------------------------

def polymer_from_cooler(
    clr,
    coordinates: str,
    insulation_scores: Optional[np.ndarray] = None,
    contact_threshold_quantile: float = 0.70,
    n_steps: int = 5000,
    dt: float = 0.005,
    friction: float = 1.0,
    temperature: float = 1.0,
    seed: int = 42,
) -> np.ndarray:
    """
    End-to-end: fetch a region from a Cooler, simulate, return coordinates.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object at the desired resolution.
    coordinates : str
        Region string, e.g. "chr12:26,000,000-28,000,000".
    insulation_scores : np.ndarray or None
        Optional per-bin insulation scores to modulate backbone stiffness.
    contact_threshold_quantile : float
        Quantile threshold for converting contacts to restraints.
    n_steps, dt, friction, temperature, seed
        Simulation parameters (see ``simulate_polymer``).

    Returns
    -------
    np.ndarray
        Bead positions, shape (N, 3).
    """
    matrix = clr.matrix(balance=True).fetch(coordinates)
    n_beads = matrix.shape[0]

    restraints = contact_matrix_to_restraints(
        matrix, contact_threshold_quantile=contact_threshold_quantile,
    )

    backbone_k = None
    if insulation_scores is not None:
        backbone_k = insulation_to_backbone_stiffness(insulation_scores)

    return simulate_polymer(
        n_beads,
        restraints,
        backbone_k=backbone_k,
        n_steps=n_steps,
        dt=dt,
        friction=friction,
        temperature=temperature,
        seed=seed,
    )

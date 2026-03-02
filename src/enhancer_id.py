"""Enhancer identification module for HiC-TAD-Library.

Pipeline
--------
1. extract_candidate_enhancers() — uses Hi-C contact structure to nominate candidates
2. one_hot_encode() / fetch_sequences() — sequence feature extraction (requires pyfaidx)
3. EnhancerCNN — lightweight dilated CNN to score candidates (requires PyTorch)
4. get_alphagenome_accessibility() — AlphaGenome proxy labels before ATAC-seq is available
5. train_enhancer_cnn() / predict_enhancers() — training and prediction

Typical workflow
----------------
    from src.tad_boundaries import (
        score_boundary_prominence, call_tad_intervals, multiscale_insulation
    )
    from src.enhancer_id import (
        extract_candidate_enhancers, fetch_sequences,
        train_enhancer_cnn, predict_enhancers,
        get_alphagenome_accessibility, EnhancerCNN,
    )

    # 1. Run existing TAD pipeline
    ins_table, windows = multiscale_insulation(clr, coords)
    boundary_df = score_boundary_prominence(ins_table, windows[2])
    tad_df = call_tad_intervals(boundary_df, clr)

    # 2. Get proxy labels from AlphaGenome
    candidates = extract_candidate_enhancers(clr, coords, tad_df, boundary_df)
    labels = get_alphagenome_accessibility(candidates, ag_client)

    # 3. Train CNN (replace labels with ATAC-seq peaks when available)
    seqs = fetch_sequences(candidates, "data/raw/mm10.fa")
    model = train_enhancer_cnn(seqs, labels)

    # 4. Predict
    enhancers = predict_enhancers(clr, coords, model, "data/raw/mm10.fa", tad_df, boundary_df)

Data dependencies
-----------------
- Cooler contact matrix (.mcool) at 5 kb resolution — same file used by other modules
- Reference genome FASTA (mm10): download to data/raw/mm10.fa, index with samtools faidx
- AlphaGenome API key in .env as ALPHA_GENOME_API_KEY
- ATAC-seq peaks (.bed) when available, to replace AlphaGenome proxy labels
"""

import os
import cooler
import numpy as np
import pandas as pd
from typing import Optional, Tuple

try:
    from src.tad_boundaries import parse_coordinates, make_view_df
except ModuleNotFoundError:
    from tad_boundaries import parse_coordinates, make_view_df  # type: ignore


# ---------------------------------------------------------------------------
# Lazy imports for optional heavy dependencies
# ---------------------------------------------------------------------------

def _import_torch() -> Tuple:
    """Lazily import PyTorch. Raises a clear error if not installed."""
    try:
        import torch
        import torch.nn as nn
        return torch, nn
    except ImportError:
        raise ImportError(
            "PyTorch is required for CNN training and prediction. "
            "Install with: conda install pytorch -c pytorch"
        )


# ---------------------------------------------------------------------------
# Stage 1: Candidate Region Extraction
# ---------------------------------------------------------------------------

def extract_candidate_enhancers(
    clr: cooler.Cooler,
    coordinates: str,
    tad_df: pd.DataFrame,
    boundary_df: pd.DataFrame,
    contact_window_bp: int = 50_000,
    boundary_exclusion_bp: int = 10_000,
) -> pd.DataFrame:
    """
    Identify candidate enhancer bins from Hi-C contact structure.

    A bin is a candidate if it:
    - Falls within a called TAD (assigned tad_id >= 0)
    - Is not within ``boundary_exclusion_bp`` of a strong boundary
    - Has an above-median local contact enrichment within the region

    The contact score per bin is the sum of balanced contacts within
    ``contact_window_bp`` on either side, excluding the two nearest diagonals
    to avoid inflating scores with self and neighbour contacts.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object at the desired resolution.
    coordinates : str
        Region string, e.g. "chr12:26,000,000-28,000,000".
    tad_df : pd.DataFrame
        Output of ``call_tad_intervals()`` — needs chrom, start, end, tad_id.
    boundary_df : pd.DataFrame
        Output of ``score_boundary_prominence()`` — needs chrom, start, end,
        boundary_class.
    contact_window_bp : int
        Window in bp for computing per-bin local contact enrichment.
    boundary_exclusion_bp : int
        Exclude bins within this distance (bp) of any strong boundary.

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, contact_score, tad_id, near_boundary.
        Rows are sorted by contact_score descending.
    """
    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    resolution = clr.binsize
    window_bins = max(1, contact_window_bp // resolution)

    matrix = clr.matrix(balance=True).fetch(coordinates)
    matrix = np.nan_to_num(matrix, nan=0.0)
    n = matrix.shape[0]

    # Per-bin contact score: off-diagonal sum within window
    contact_scores = np.zeros(n, dtype=np.float32)
    for i in range(n):
        lo = max(0, i - window_bins)
        hi = min(n, i + window_bins + 1)
        row = matrix[i, lo:hi].copy()
        diag_lo = i - lo
        # Zero out ±1 diagonal to avoid self/neighbour contact inflation
        for d in (-1, 0, 1):
            j = diag_lo + d
            if 0 <= j < len(row):
                row[j] = 0.0
        contact_scores[i] = float(row.sum())

    bins = clr.bins().fetch(coordinates)
    candidates = bins[['chrom', 'start', 'end']].iloc[:n].copy().reset_index(drop=True)
    candidates['contact_score'] = contact_scores

    # Assign TAD membership
    candidates['tad_id'] = -1
    for _, tad_row in tad_df[tad_df['chrom'] == chrom].iterrows():
        mask = (
            (candidates['start'] >= tad_row['start']) &
            (candidates['end'] <= tad_row['end'])
        )
        candidates.loc[mask, 'tad_id'] = int(tad_row['tad_id'])

    # Mark bins near strong boundaries
    candidates['near_boundary'] = False
    strong_bounds = boundary_df[
        (boundary_df['chrom'] == chrom) &
        (boundary_df['boundary_class'] == 'strong')
    ]
    bin_mids = candidates['start'].values + resolution // 2
    for _, bnd in strong_bounds.iterrows():
        bnd_mid = int((bnd['start'] + bnd['end']) // 2)
        candidates.loc[np.abs(bin_mids - bnd_mid) < boundary_exclusion_bp, 'near_boundary'] = True

    # Filter: inside a TAD and not adjacent to a strong boundary
    filtered = candidates[
        (candidates['tad_id'] >= 0) & (~candidates['near_boundary'])
    ].copy()

    # Keep above-median contact score to reduce noise
    if len(filtered) > 0:
        median_score = filtered['contact_score'].median()
        filtered = filtered[filtered['contact_score'] >= median_score].copy()

    return (
        filtered
        .sort_values('contact_score', ascending=False)
        .reset_index(drop=True)
    )


# ---------------------------------------------------------------------------
# Stage 2: Sequence Feature Extraction
# ---------------------------------------------------------------------------

_ENCODE_MAP: dict = {
    'A': 0, 'C': 1, 'G': 2, 'T': 3,
    'a': 0, 'c': 1, 'g': 2, 't': 3,
}


def one_hot_encode(sequence: str) -> np.ndarray:
    """
    One-hot encode a DNA sequence.

    Parameters
    ----------
    sequence : str
        DNA string (ACGT, any case). Ambiguous bases (N, etc.) → all zeros.

    Returns
    -------
    np.ndarray
        Shape (4, len(sequence)), dtype float32. Row order: A, C, G, T.
    """
    L = len(sequence)
    encoded = np.zeros((4, L), dtype=np.float32)
    for j, base in enumerate(sequence):
        idx = _ENCODE_MAP.get(base)
        if idx is not None:
            encoded[idx, j] = 1.0
    return encoded


def fetch_sequences(
    candidate_df: pd.DataFrame,
    fasta_path: str,
    window_bp: int = 500,
) -> np.ndarray:
    """
    Fetch and one-hot encode DNA sequences centered on each candidate.

    Requires ``pyfaidx`` and an indexed FASTA file (.fai).
    Index with: ``samtools faidx data/raw/mm10.fa``

    Parameters
    ----------
    candidate_df : pd.DataFrame
        Requires chrom, start, end columns (e.g. output of
        ``extract_candidate_enhancers()``).
    fasta_path : str
        Path to reference genome FASTA file.
    window_bp : int
        Total sequence window in bp, centered on each candidate midpoint.

    Returns
    -------
    np.ndarray
        Shape (n_candidates, 4, window_bp), dtype float32.
    """
    try:
        from pyfaidx import Fasta
    except ImportError:
        raise ImportError(
            "pyfaidx is required for sequence fetching. "
            "Install with: pip install pyfaidx"
        )

    genome = Fasta(fasta_path)
    half = window_bp // 2
    results = []

    for _, row in candidate_df.iterrows():
        mid = int((row['start'] + row['end']) // 2)
        seq_start = max(0, mid - half)
        seq_end = seq_start + window_bp
        try:
            seq = str(genome[row['chrom']][seq_start:seq_end])
        except (KeyError, ValueError):
            seq = 'N' * window_bp

        # Pad to exact window_bp if near chromosome end
        if len(seq) < window_bp:
            seq = seq + 'N' * (window_bp - len(seq))

        results.append(one_hot_encode(seq[:window_bp]))

    return np.stack(results, axis=0)  # (n_candidates, 4, window_bp)


# ---------------------------------------------------------------------------
# Stage 3: CNN Architecture
# ---------------------------------------------------------------------------

class EnhancerCNN:
    """
    Lightweight dilated CNN for predicting enhancer activity from DNA sequence.

    Requires PyTorch (lazily imported at instantiation). The CNN is deliberately
    shallow to train on the modest number of candidates available from a single
    genomic region, avoiding overfitting.

    Architecture
    ------------
    Input (batch, 4, sequence_length)
        Conv1d(4→n_filters, k=kernel_size, pad='same') → BN → ReLU
        Conv1d(n_filters→n_filters, k=kernel_size, dilation=2, pad='same') → BN → ReLU
        Conv1d(n_filters→n_filters, k=kernel_size, dilation=4, pad='same') → BN → ReLU
        GlobalAvgPool1d → Linear(n_filters→1) → Sigmoid
    Output (batch, 1)  enhancer score ∈ [0, 1]

    The dilated convolutions capture motifs at three scales:
    - dilation=1: local TF binding motifs (~8 bp)
    - dilation=2: short-range motif grammar (~16 bp)
    - dilation=4: composite regulatory elements (~32 bp)

    Usage
    -----
    model = train_enhancer_cnn(sequences, labels)
    enhancers = predict_enhancers(clr, coords, model, fasta_path, tad_df, boundary_df)
    """

    def __init__(self, n_filters: int = 64, kernel_size: int = 8):
        torch, nn = _import_torch()
        self._torch = torch

        class _CNN(nn.Module):
            def __init__(self):
                super().__init__()
                self.conv1 = nn.Sequential(
                    nn.Conv1d(4, n_filters, kernel_size, padding='same'),
                    nn.BatchNorm1d(n_filters),
                    nn.ReLU(),
                )
                self.conv2 = nn.Sequential(
                    nn.Conv1d(n_filters, n_filters, kernel_size,
                              dilation=2, padding='same'),
                    nn.BatchNorm1d(n_filters),
                    nn.ReLU(),
                )
                self.conv3 = nn.Sequential(
                    nn.Conv1d(n_filters, n_filters, kernel_size,
                              dilation=4, padding='same'),
                    nn.BatchNorm1d(n_filters),
                    nn.ReLU(),
                )
                self.head = nn.Sequential(
                    nn.Linear(n_filters, 1),
                    nn.Sigmoid(),
                )

            def forward(self, x):
                x = self.conv1(x)
                x = self.conv2(x)
                x = self.conv3(x)
                x = x.mean(dim=-1)  # Global average pool over sequence length
                return self.head(x)

        self._model = _CNN()

    @property
    def model(self):
        return self._model

    def save(self, path: str) -> None:
        """Save model weights to a .pt file."""
        self._torch.save(self._model.state_dict(), path)

    def load(self, path: str) -> None:
        """Load model weights from a .pt file."""
        self._model.load_state_dict(
            self._torch.load(path, map_location='cpu', weights_only=True)
        )
        self._model.eval()


# ---------------------------------------------------------------------------
# Stage 4: AlphaGenome Accessibility Labels
# ---------------------------------------------------------------------------

def get_alphagenome_accessibility(
    candidate_df: pd.DataFrame,
    ag_client,
    track_name: str = 'ATAC',
) -> np.ndarray:
    """
    Query AlphaGenome for predicted accessibility at each candidate locus.

    Returns proxy training labels normalized to [0, 1]. Use these before
    real ATAC-seq data is available; swap in ATAC-seq peak scores once
    the experiment is run.

    Parameters
    ----------
    candidate_df : pd.DataFrame
        Requires chrom, start, end columns.
    ag_client
        Initialized AlphaGenome client (from src/alphagenome/).
    track_name : str
        Accessibility track to query ('ATAC', 'DNase', or histone track).

    Returns
    -------
    np.ndarray
        Shape (n_candidates,), float32. Normalized accessibility scores.
    """
    scores = np.zeros(len(candidate_df), dtype=np.float32)

    for idx, (_, row) in enumerate(candidate_df.iterrows()):
        mid = int((row['start'] + row['end']) // 2)
        try:
            result = ag_client.predict(
                chrom=row['chrom'],
                pos=mid,
                tracks=[track_name],
            )
            if result is not None and track_name in result:
                track_vals = np.asarray(result[track_name], dtype=np.float32)
                center = len(track_vals) // 2
                # Average over ±2 bins around the center (AlphaGenome is 128bp/bin)
                scores[idx] = float(track_vals[max(0, center - 2):center + 3].mean())
        except Exception:
            scores[idx] = 0.0

    # Normalize to [0, 1]
    lo, hi = scores.min(), scores.max()
    if hi > lo:
        scores = (scores - lo) / (hi - lo)

    return scores


# ---------------------------------------------------------------------------
# Stage 5: Training and Prediction
# ---------------------------------------------------------------------------

def train_enhancer_cnn(
    sequences: np.ndarray,
    labels: np.ndarray,
    epochs: int = 20,
    lr: float = 1e-3,
    batch_size: int = 32,
    save_path: Optional[str] = None,
    verbose: bool = True,
) -> 'EnhancerCNN':
    """
    Train the EnhancerCNN on (sequence, accessibility) pairs.

    Parameters
    ----------
    sequences : np.ndarray
        Shape (N, 4, window_bp). Output of ``fetch_sequences()``.
    labels : np.ndarray
        Shape (N,). Accessibility scores in [0, 1].
        Source: ``get_alphagenome_accessibility()`` or ATAC-seq peak scores.
    epochs : int
        Number of training epochs.
    lr : float
        Adam optimizer learning rate.
    batch_size : int
        Mini-batch size.
    save_path : str or None
        If provided, saves model weights here after training.
        Example: "data/processed/enhancer_cnn.pt"
    verbose : bool
        Print loss every 5 epochs.

    Returns
    -------
    EnhancerCNN
        Trained model in eval mode.
    """
    torch, nn = _import_torch()

    model = EnhancerCNN()
    optimizer = torch.optim.Adam(model.model.parameters(), lr=lr)
    criterion = nn.MSELoss()

    X = torch.from_numpy(sequences).float()
    y = torch.from_numpy(labels.astype(np.float32)).unsqueeze(1)

    model.model.train()
    n = len(X)

    for epoch in range(epochs):
        perm = torch.randperm(n)
        epoch_loss = 0.0
        n_batches = 0
        for i in range(0, n, batch_size):
            idx_b = perm[i:i + batch_size]
            xb, yb = X[idx_b], y[idx_b]
            optimizer.zero_grad()
            loss = criterion(model.model(xb), yb)
            loss.backward()
            optimizer.step()
            epoch_loss += loss.item()
            n_batches += 1

        if verbose and (epoch + 1) % 5 == 0:
            print(f"  Epoch {epoch + 1}/{epochs}  loss={epoch_loss / max(n_batches, 1):.4f}")

    model.model.eval()

    if save_path is not None:
        dir_path = os.path.dirname(save_path)
        if dir_path:
            os.makedirs(dir_path, exist_ok=True)
        model.save(save_path)

    return model


def predict_enhancers(
    clr: cooler.Cooler,
    coordinates: str,
    model: 'EnhancerCNN',
    fasta_path: str,
    tad_df: pd.DataFrame,
    boundary_df: pd.DataFrame,
    window_bp: int = 500,
    threshold: float = 0.5,
    contact_window_bp: int = 50_000,
    boundary_exclusion_bp: int = 10_000,
) -> pd.DataFrame:
    """
    End-to-end enhancer prediction for a genomic region.

    Chains: candidate extraction → sequence fetch → CNN scoring.

    Parameters
    ----------
    clr : cooler.Cooler
        Cooler object at the desired resolution.
    coordinates : str
        Region string, e.g. "chr12:26,000,000-28,000,000".
    model : EnhancerCNN
        Trained CNN model (from ``train_enhancer_cnn()``).
    fasta_path : str
        Path to reference genome FASTA (must be indexed with .fai).
    tad_df : pd.DataFrame
        Output of ``call_tad_intervals()``.
    boundary_df : pd.DataFrame
        Output of ``score_boundary_prominence()``.
    window_bp : int
        Sequence window in bp for CNN input.
    threshold : float
        enhancer_score threshold for the ``is_enhancer`` boolean flag.
    contact_window_bp : int
        Window for per-bin contact score calculation.
    boundary_exclusion_bp : int
        Exclusion zone around strong boundaries in bp.

    Returns
    -------
    pd.DataFrame
        Columns: chrom, start, end, enhancer_score, is_enhancer,
        contact_score, tad_id, candidate_rank.
        Sorted by enhancer_score descending.
    """
    torch, _ = _import_torch()

    candidates = extract_candidate_enhancers(
        clr, coordinates, tad_df, boundary_df,
        contact_window_bp=contact_window_bp,
        boundary_exclusion_bp=boundary_exclusion_bp,
    )

    if len(candidates) == 0:
        return pd.DataFrame(columns=[
            'chrom', 'start', 'end', 'enhancer_score', 'is_enhancer',
            'contact_score', 'tad_id', 'candidate_rank',
        ])

    sequences = fetch_sequences(candidates, fasta_path, window_bp=window_bp)

    model.model.eval()
    with torch.no_grad():
        X = torch.from_numpy(sequences).float()
        scores = model.model(X).squeeze(1).cpu().numpy().astype(np.float32)

    candidates = candidates.copy()
    candidates['enhancer_score'] = scores
    candidates['is_enhancer'] = scores >= threshold
    candidates['candidate_rank'] = (
        candidates['enhancer_score'].rank(ascending=False, method='min').astype(int)
    )

    return (
        candidates[[
            'chrom', 'start', 'end', 'enhancer_score', 'is_enhancer',
            'contact_score', 'tad_id', 'candidate_rank',
        ]]
        .sort_values('enhancer_score', ascending=False)
        .reset_index(drop=True)
    )

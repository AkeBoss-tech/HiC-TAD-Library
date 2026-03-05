"""
Synthetic enhancer design module for HiC-TAD-Library.

Two complementary design strategies:

1. Gradient-based sequence optimization
   Backpropagate through the EnhancerCNN into a soft (Gumbel-softmax) one-hot
   representation of the sequence. Starting from random logits, iteratively
   update toward sequences that maximise the CNN score.

2. Motif-insertion combinatorial search
   Take a seed sequence (e.g. the top candidate from predict_enhancers),
   systematically try inserting known neural TF motifs at every position,
   and greedily accept insertions that increase CNN score.

Both strategies require a trained EnhancerCNN (src/enhancer_id.py).
PyTorch is imported lazily at call time.

Typical workflow
----------------
    from src.enhancer_id import EnhancerCNN
    from src.enhancer_design import gradient_design, motif_insert_design, score_sequences

    cnn = EnhancerCNN()
    cnn.load("data/processed/enhancer_cnn.pt")

    # Generate 10 synthetic enhancers by gradient optimisation
    seqs, scores = gradient_design(cnn, n_sequences=10, steps=500)

    # Or refine a real candidate sequence via motif insertion
    improved, score = motif_insert_design(cnn, seed_sequence)

    # AlphaGenome in-silico validation (needs real locus context)
    results = score_sequences_alphagenome(ag_model, locus_interval, locus_seq, designs)
"""

import numpy as np
from typing import List, Tuple, Optional, Dict

# Neural TF motifs relevant for Sox11 / Mir9-2 regulatory context
# Each value is one or more consensus sequences for that factor.
# Sources: JASPAR2024 vertebrate motif database (top-information-content consensus).
NEURAL_MOTIFS: Dict[str, List[str]] = {
    "Sox11":    ["AACAAAG", "AACAATG"],   # Sox family HMG-box consensus
    "Sox2":     ["CTTTGTT", "CTTTGTC"],
    "NeuroD1":  ["CAGATGG"],              # bHLH E-box variant
    "Atoh1":    ["CAGCTGT"],              # bHLH E-box variant
    "AP1":      ["TGAGTCA", "TGACTCA"],   # Fos/Jun composite
    "KLF":      ["CCCCGCCC"],             # KLF4/KLF9 GC-box
    "CTCF":     ["CCGCGAGGCGCAGG"],       # CTCF 19-mer (include to test boundary effect)
    "Ebox":     ["CACGTG", "CATGTG"],     # Generic E-box
    "SP1":      ["GGGCGG"],               # SP1/SP3 GC-rich
    "NRSF":     ["TTCAGCACCACGGACAGCGCC"], # REST/NRSF repressor — use as negative control
}


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _import_torch():
    try:
        import torch
        import torch.nn as nn
        import torch.nn.functional as F
        return torch, nn, F
    except ImportError:
        raise ImportError(
            "PyTorch is required for sequence design. "
            "Install with: conda install pytorch -c pytorch"
        )


_BASES = ['A', 'C', 'G', 'T']
_BASE_IDX = {b: i for i, b in enumerate(_BASES)}


def _seq_to_tensor(seq: str):
    """Convert ACGT string → (1, 4, L) float32 tensor."""
    torch, _, _ = _import_torch()
    L = len(seq)
    t = torch.zeros(1, 4, L, dtype=torch.float32)
    for j, b in enumerate(seq.upper()):
        idx = _BASE_IDX.get(b)
        if idx is not None:
            t[0, idx, j] = 1.0
    return t


def _tensor_to_seq(t) -> str:
    """Convert (4, L) or (1, 4, L) float tensor → ACGT string (argmax decode)."""
    if t.dim() == 3:
        t = t.squeeze(0)
    indices = t.argmax(dim=0).tolist()
    return ''.join(_BASES[i] for i in indices)


def gc_content(seq: str) -> float:
    """Fraction of G+C bases."""
    seq = seq.upper()
    if not seq:
        return 0.0
    return (seq.count('G') + seq.count('C')) / len(seq)


# ---------------------------------------------------------------------------
# Strategy 1: Gradient-based sequence optimisation
# ---------------------------------------------------------------------------

def gradient_design(
    cnn_model,
    n_sequences: int = 10,
    length: int = 500,
    steps: int = 600,
    lr: float = 0.05,
    gc_penalty: float = 0.5,
    init_gc: float = 0.45,
    verbose: bool = True,
) -> Tuple[List[str], List[float]]:
    """
    Design synthetic enhancer sequences by gradient-based optimisation.

    Uses temperature-annealed Gumbel-softmax to maintain a differentiable
    approximation to discrete sequence identity. Includes a GC-content
    penalty to discourage degenerate GC-rich or AT-rich solutions.

    Parameters
    ----------
    cnn_model : EnhancerCNN
        A trained EnhancerCNN instance (from src.enhancer_id).
    n_sequences : int
        Number of independent designs to generate.
    length : int
        Sequence length in bp (should match CNN training window, default 500).
    steps : int
        Gradient optimisation steps per sequence.
    lr : float
        Adam learning rate for logit updates.
    gc_penalty : float
        Weight of L2 penalty toward init_gc GC content. 0 = no penalty.
    init_gc : float
        Target GC fraction for penalty (0.45 = typical mammalian enhancer).
    verbose : bool
        Print progress every 100 steps.

    Returns
    -------
    sequences : List[str]
        Designed ACGT sequences.
    scores : List[float]
        CNN scores ∈ [0, 1] for each designed sequence.
    """
    torch, nn, F = _import_torch()

    model = cnn_model.model
    model.eval()
    for p in model.parameters():
        p.requires_grad_(False)

    sequences, scores = [], []

    for i in range(n_sequences):
        # Initialise logits with slight GC bias matching mammalian genome
        # A/T logits slightly lower → tends toward ~45% GC at softmax
        logits = torch.zeros(1, 4, length)
        logits[0, 1, :] = 0.1  # C
        logits[0, 2, :] = 0.1  # G
        logits = logits + 0.3 * torch.randn_like(logits)
        logits.requires_grad_(True)

        optimizer = torch.optim.Adam([logits], lr=lr)
        # Temperature schedule: start hot (1.0), anneal to cold (0.1)
        temps = np.linspace(1.0, 0.1, steps)

        best_score = -1.0
        best_logits = logits.detach().clone()

        for step in range(steps):
            optimizer.zero_grad()
            tau = float(temps[step])
            soft_seq = F.gumbel_softmax(logits, tau=tau, hard=False, dim=1)

            score = model(soft_seq)  # (1, 1)

            # GC penalty: push GC fraction toward init_gc
            gc_frac = (soft_seq[0, 1, :] + soft_seq[0, 2, :]).mean()
            penalty = gc_penalty * (gc_frac - init_gc) ** 2

            loss = -score.mean() + penalty
            loss.backward()
            optimizer.step()

            s = float(score.item())
            if s > best_score:
                best_score = s
                best_logits = logits.detach().clone()

        # Decode best logits
        with torch.no_grad():
            seq = _tensor_to_seq(best_logits.squeeze(0))
            # Compute final hard score
            hard_t = _seq_to_tensor(seq)
            final_score = float(model(hard_t).item())

        sequences.append(seq)
        scores.append(final_score)

        if verbose:
            print(f"  Design {i+1}/{n_sequences}: score={final_score:.4f}  GC={gc_content(seq):.2%}")

    return sequences, scores


# ---------------------------------------------------------------------------
# Strategy 2: Motif-insertion combinatorial search
# ---------------------------------------------------------------------------

def motif_insert_design(
    cnn_model,
    seed_sequence: str,
    motif_set: Optional[Dict[str, List[str]]] = None,
    n_rounds: int = 3,
    verbose: bool = True,
) -> Tuple[str, float, List[dict]]:
    """
    Improve a seed sequence by greedily inserting neural TF motifs.

    For each motif, tries inserting at every valid position. Accepts the
    best single insertion that increases CNN score, then repeats for
    n_rounds. This is a greedy hill-climb in the space of motif-augmented
    sequences.

    Parameters
    ----------
    cnn_model : EnhancerCNN
        A trained EnhancerCNN instance.
    seed_sequence : str
        Starting ACGT sequence (typically the best real candidate).
    motif_set : dict, optional
        {factor_name: [consensus, ...]}. Defaults to NEURAL_MOTIFS.
    n_rounds : int
        Number of greedy insertion rounds.
    verbose : bool
        Print each accepted insertion.

    Returns
    -------
    best_seq : str
        Final sequence after all rounds.
    best_score : float
        CNN score of final sequence.
    history : List[dict]
        Record of each accepted insertion: {round, factor, motif, position, score}.
    """
    torch, _, _ = _import_torch()

    if motif_set is None:
        motif_set = NEURAL_MOTIFS

    model = cnn_model.model
    model.eval()

    def _score(seq: str) -> float:
        with torch.no_grad():
            return float(model(_seq_to_tensor(seq)).item())

    current_seq = seed_sequence.upper()
    current_score = _score(current_seq)
    history = []
    L = len(current_seq)

    if verbose:
        print(f"  Seed score: {current_score:.4f}  (GC={gc_content(current_seq):.2%})")

    for rnd in range(n_rounds):
        best_delta = 0.0
        best_variant = None
        best_meta = None

        for factor, motifs in motif_set.items():
            for motif in motifs:
                m_len = len(motif)
                if m_len >= L:
                    continue
                # Try each valid insertion position
                step = max(1, m_len // 2)
                for pos in range(0, L - m_len, step):
                    variant = current_seq[:pos] + motif + current_seq[pos + m_len:]
                    variant = variant[:L]  # keep length fixed
                    s = _score(variant)
                    delta = s - current_score
                    if delta > best_delta:
                        best_delta = delta
                        best_variant = variant
                        best_meta = {"round": rnd + 1, "factor": factor,
                                     "motif": motif, "position": pos, "score": s}

        if best_variant is not None:
            current_seq = best_variant
            current_score = best_meta["score"]
            history.append(best_meta)
            if verbose:
                print(f"  Round {rnd+1}: insert {best_meta['factor']} "
                      f"@ pos {best_meta['position']} → score {current_score:.4f} "
                      f"(+{best_delta:.4f})")
        else:
            if verbose:
                print(f"  Round {rnd+1}: no improvement found — stopping early.")
            break

    return current_seq, current_score, history


# ---------------------------------------------------------------------------
# Scoring utilities
# ---------------------------------------------------------------------------

def score_sequences(cnn_model, sequences: List[str]) -> np.ndarray:
    """
    Score a list of ACGT sequences with the trained CNN.

    Parameters
    ----------
    cnn_model : EnhancerCNN
    sequences : List[str]
        Each sequence should be length == CNN training window (500 bp).

    Returns
    -------
    np.ndarray
        Shape (n_sequences,), CNN scores ∈ [0, 1].
    """
    torch, _, _ = _import_torch()
    model = cnn_model.model
    model.eval()

    scores = []
    with torch.no_grad():
        for seq in sequences:
            t = _seq_to_tensor(seq)
            scores.append(float(model(t).item()))
    return np.array(scores, dtype=np.float32)


def summarise_designs(
    sequences: List[str],
    scores: np.ndarray,
    labels: Optional[List[str]] = None,
) -> "pd.DataFrame":
    """
    Build a summary DataFrame for designed sequences.

    Columns: label, length, gc_content, cnn_score, sequence.
    """
    import pandas as pd

    if labels is None:
        labels = [f"design_{i+1}" for i in range(len(sequences))]

    rows = []
    for lbl, seq, sc in zip(labels, sequences, scores):
        rows.append({
            "label":      lbl,
            "length":     len(seq),
            "gc_content": round(gc_content(seq), 4),
            "cnn_score":  round(float(sc), 4),
            "sequence":   seq,
        })
    return pd.DataFrame(rows).sort_values("cnn_score", ascending=False).reset_index(drop=True)


# ---------------------------------------------------------------------------
# AlphaGenome in-silico validation
# ---------------------------------------------------------------------------

def validate_designs_alphagenome(
    ag_model,
    locus_sequence: str,
    locus_interval,
    candidate_start: int,
    designs: List[str],
    ontology_terms: List[str],
) -> "pd.DataFrame":
    """
    Validate designed sequences by substituting them into a real locus
    context and querying AlphaGenome for predicted ATAC change.

    Computes delta-ATAC = predicted_ATAC(design) − predicted_ATAC(wildtype)
    at the candidate bin position.

    Parameters
    ----------
    ag_model : DnaClient
        Initialised AlphaGenome client.
    locus_sequence : str
        Full ~1 Mb (2^20 bp) wildtype sequence for the locus.
    locus_interval : alphagenome.data.genome.Interval
        Genomic interval corresponding to locus_sequence.
    candidate_start : int
        Start position (within locus_sequence) of the 500 bp window to replace.
    designs : List[str]
        Designed sequences (each 500 bp).
    ontology_terms : List[str]
        ATAC ontology terms to query (e.g. ATAC_ONTOLOGIES from train_enhancers).

    Returns
    -------
    pd.DataFrame
        Columns: design_index, delta_atac_mean, delta_atac_max, sequence.

    Notes
    -----
    AlphaGenome predict_sequence() requires exactly 2^20 bp (1,048,576 bp).
    locus_sequence must match this length. The designed 500 bp window is
    substituted at candidate_start:candidate_start+500.
    """
    import pandas as pd
    from alphagenome.models import dna_model, dna_output

    REQUIRED_LEN = 2 ** 20
    if len(locus_sequence) != REQUIRED_LEN:
        raise ValueError(
            f"locus_sequence must be exactly {REQUIRED_LEN} bp, got {len(locus_sequence)} bp."
        )
    window = 500
    if candidate_start + window > REQUIRED_LEN:
        raise ValueError("candidate_start + 500 exceeds locus_sequence length.")

    # Wildtype prediction
    wt_output = ag_model.predict_sequence(
        locus_sequence,
        organism=dna_model.Organism.MUS_MUSCULUS,
        requested_outputs=[dna_output.OutputType.ATAC],
        ontology_terms=ontology_terms,
        interval=locus_interval,
    )
    wt_atac = wt_output.atac.values.mean(axis=1)  # (n_positions,) averaged over tracks

    rows = []
    for i, design_seq in enumerate(designs):
        # Substitute designed sequence into locus
        mut_seq = (
            locus_sequence[:candidate_start]
            + design_seq[:window]
            + locus_sequence[candidate_start + window:]
        )
        mut_output = ag_model.predict_sequence(
            mut_seq,
            organism=dna_model.Organism.MUS_MUSCULUS,
            requested_outputs=[dna_output.OutputType.ATAC],
            ontology_terms=ontology_terms,
            interval=locus_interval,
        )
        mut_atac = mut_output.atac.values.mean(axis=1)

        delta = mut_atac - wt_atac
        # Summarise delta around the substituted region (±5 kb)
        n_pos = len(delta)
        region_start_idx = max(0, (candidate_start // 128) - 40)
        region_end_idx   = min(n_pos, (candidate_start // 128) + 40)
        local_delta = delta[region_start_idx:region_end_idx]

        rows.append({
            "design_index":    i + 1,
            "delta_atac_mean": float(local_delta.mean()),
            "delta_atac_max":  float(local_delta.max()),
            "sequence":        design_seq[:window],
        })

    return pd.DataFrame(rows).sort_values("delta_atac_mean", ascending=False).reset_index(drop=True)

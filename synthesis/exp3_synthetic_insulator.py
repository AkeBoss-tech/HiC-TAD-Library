#!/usr/bin/env python3
"""
Experiment 3 — Synthetic stronger insulator.

Design and score CTCF insertions that beat wild-type insulation at Unc5b.

5 insertion variants:
  1pair_at_gap, 2pair_at_gap, 3pair_at_gap  (gap = lowest-signal position in cluster)
  1pair_5kb_left, 1pair_5kb_right           (flanking positions)

Output:
    media/synthesis_exp3_synthetic_insulator.png
"""

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
from dotenv import load_dotenv

# ── Constants ─────────────────────────────────────────────────────────────────
CHROM     = "chr10"
WIN_START = 60_240_000
WIN_END   = 61_288_576
CELL_TYPE = "EFO:0004038"
WT_CACHE  = "data/processed/synthesis_unc5b_wt_cache.npz"
OUT_PNG   = "media/synthesis_exp3_synthetic_insulator.png"

CTCF_FWD  = "CCGCGAGGCGCAGG"
CTCF_REV  = "CCTGCGCCTCGCGG"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_model():
    load_dotenv()
    from alphagenome.models import dna_client
    api_key = os.getenv("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise EnvironmentError("ALPHA_GENOME_API_KEY not set in .env")
    return dna_client.create(api_key), dna_client


def compute_insulation(mat, resolution, window_bp=50_000):
    n   = mat.shape[0]
    win = max(2, int(window_bp / resolution))
    ins = np.full(n, np.nan)
    for i in range(win, n - win):
        ins[i] = mat[i - win:i, i:i + win].mean()
    return ins


def load_wt_cache():
    if not os.path.exists(WT_CACHE):
        return None
    try:
        data = np.load(WT_CACHE)
        return dict(
            contact_mat=data["contact_mat"],
            ctcf_signal=data["ctcf_signal"],
            resolution=float(data["resolution"]),
            n_bins=int(data["n_bins"]),
        )
    except Exception as exc:
        warnings.warn(f"WT cache corrupt ({exc}); will recompute.")
        return None


def predict_wildtype(model, dna_client):
    from alphagenome.data import genome
    interval = genome.Interval(chromosome=CHROM, start=WIN_START, end=WIN_END)
    organism = dna_client.Organism.MUS_MUSCULUS
    print(f"Predicting WT: {CHROM}:{WIN_START:,}–{WIN_END:,} ...")
    wt_out = model.predict_interval(
        interval=interval, organism=organism,
        requested_outputs={
            dna_client.OutputType.CONTACT_MAPS,
            dna_client.OutputType.CHIP_TF,
        },
        ontology_terms=[CELL_TYPE],
    )
    contact_mat = wt_out.contact_maps.values[:, :, 0]
    n_bins      = contact_mat.shape[0]
    resolution  = (WIN_END - WIN_START) / n_bins
    try:
        tf_meta     = wt_out.chip_tf.metadata
        ctcf_mask   = tf_meta["transcription_factor"] == "CTCF"
        ctcf_signal = (
            wt_out.chip_tf.values[:, ctcf_mask].mean(axis=1)
            if ctcf_mask.sum() > 0
            else wt_out.chip_tf.values.mean(axis=1)
        )
    except Exception:
        ctcf_signal = wt_out.chip_tf.values.mean(axis=1)
    return dict(contact_mat=contact_mat, ctcf_signal=ctcf_signal,
                resolution=resolution, n_bins=n_bins)


def build_ctcf_pair_seq(n_pairs, spacer=10):
    """Build convergent CTCF pair sequence (FWD + spacer + REV, repeated)."""
    pair = CTCF_FWD + "N" * spacer + CTCF_REV
    return ("N" * spacer).join([pair] * n_pairs)


def find_cluster_bounds(ctcf_signal):
    """Return (cluster_start, cluster_end) as indices of first/last above-median bin."""
    median = np.nanmedian(ctcf_signal)
    above  = np.where(ctcf_signal > median)[0]
    if len(above) < 2:
        n = len(ctcf_signal)
        return n // 4, 3 * n // 4
    return int(above[0]), int(above[-1])


def find_gap(ctcf_signal, cluster_start, cluster_end):
    """Lowest-signal position between first and last CTCF peak."""
    region = ctcf_signal[cluster_start:cluster_end]
    gap_offset = int(np.argmin(region))
    return cluster_start + gap_offset


def validate_position(pos, label):
    if not (WIN_START <= pos < WIN_END):
        raise ValueError(
            f"Variant position {pos:,} out of window "
            f"{WIN_START:,}–{WIN_END:,} for '{label}'"
        )


def predict_variant_safe(model, dna_client, interval, variant, label):
    try:
        organism = dna_client.Organism.MUS_MUSCULUS
        var_out  = model.predict_variant(
            interval=interval, variant=variant, organism=organism,
            requested_outputs={dna_client.OutputType.CONTACT_MAPS},
            ontology_terms=[CELL_TYPE],
        )
        return var_out.alternate.contact_maps.values[:, :, 0]
    except Exception as exc:
        warnings.warn(f"Variant '{label}' failed: {exc}")
        return None


def boundary_bin(resolution):
    return int((WIN_END - WIN_START) / 2 / resolution)


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    model, dna_client = load_model()
    from alphagenome.data import genome

    wt_data = load_wt_cache()
    if wt_data is None:
        wt_data = predict_wildtype(model, dna_client)

    contact_mat = wt_data["contact_mat"]
    ctcf_signal = wt_data["ctcf_signal"]
    resolution  = wt_data["resolution"]
    n_bins      = wt_data["n_bins"]

    chip_tf_res  = (WIN_END - WIN_START) / len(ctcf_signal)
    cluster_s, cluster_e = find_cluster_bounds(ctcf_signal)
    gap_idx      = find_gap(ctcf_signal, cluster_s, cluster_e)
    gap_genomic  = WIN_START + int(gap_idx * chip_tf_res)
    print(f"Gap position: bin {gap_idx} → {gap_genomic:,}")

    interval = genome.Interval(chromosome=CHROM, start=WIN_START, end=WIN_END)

    # Build 5 insertion variants
    designs = {}
    for n_pairs in [1, 2, 3]:
        motif = build_ctcf_pair_seq(n_pairs)
        label = f"{n_pairs}pair_at_gap"
        validate_position(gap_genomic, label)
        designs[label] = genome.Variant(
            chromosome=CHROM, position=gap_genomic - 1,
            reference_bases="N" * len(motif), alternate_bases=motif,
            name=label,
        )

    motif_1 = build_ctcf_pair_seq(1)
    for shift, suffix in [(-5000, "5kb_left"), (+5000, "5kb_right")]:
        pos   = gap_genomic + shift
        label = f"1pair_{suffix}"
        # Clamp to window
        pos = max(WIN_START, min(WIN_END - len(motif_1), pos))
        validate_position(pos, label)
        designs[label] = genome.Variant(
            chromosome=CHROM, position=pos - 1,
            reference_bases="N" * len(motif_1), alternate_bases=motif_1,
            name=label,
        )

    # Predict all designs
    ins_wt   = compute_insulation(contact_mat, resolution)
    bbin     = boundary_bin(resolution)
    results  = {}
    mut_mats = {}

    for label, variant in designs.items():
        print(f"Predicting: {label} ...")
        mat = predict_variant_safe(model, dna_client, interval, variant, label)
        if mat is None:
            mat = contact_mat.copy()
        mut_mats[label] = mat
        ins_mut = compute_insulation(mat, resolution)
        delta   = float(ins_mut[bbin] - ins_wt[bbin])
        results[label] = dict(ins=ins_mut, delta=delta)
        print(f"  Δ insulation at boundary: {delta:+.5f}")

    # Rank by delta descending (positive = stronger insulation)
    ranked = sorted(results.items(), key=lambda x: x[1]["delta"], reverse=True)
    print("\nRanking (best to worst):")
    for i, (lbl, r) in enumerate(ranked):
        print(f"  {i+1}. {lbl}: Δ={r['delta']:+.5f}")

    best_label = ranked[0][0]

    # ── Figure: 2-panel ───────────────────────────────────────────────────────
    xbins_mb = np.linspace(WIN_START / 1e6, WIN_END / 1e6, n_bins)
    vmax     = np.nanpercentile(contact_mat[contact_mat > 0], 99)
    extent   = [WIN_START / 1e6, WIN_END / 1e6, WIN_START / 1e6, WIN_END / 1e6]

    colors = plt.cm.tab10(np.linspace(0, 0.9, len(designs)))

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Panel 1: insulation score lines
    ax1.plot(xbins_mb, ins_wt, color="black", lw=2, linestyle="--", label="WT", zorder=10)
    for (label, r), c in zip(results.items(), colors):
        ax1.plot(xbins_mb, r["ins"], color=c, lw=1.5,
                 label=f"{label} (Δ={r['delta']:+.4f})")
    bound_mb = WIN_START / 1e6 + (WIN_END - WIN_START) / 2 / 1e6
    ax1.axvline(bound_mb, color="red", lw=1.5, linestyle="-", alpha=0.8, label="Boundary")
    ax1.set_xlabel("Genomic position (Mb)")
    ax1.set_ylabel("Insulation score")
    ax1.set_title("Insulation Score — WT vs Synthetic Designs", fontweight="bold")
    ax1.legend(fontsize=7, loc="lower right")

    # Panel 2: side-by-side contact maps (WT left, best design right)
    n_half = n_bins // 2
    combined = np.zeros((n_bins, n_bins))
    combined[:, :n_half] = contact_mat[:, :n_half]
    combined[:, n_half:] = mut_mats[best_label][:, n_half:]
    ax2.imshow(combined, cmap="Reds", aspect="auto", origin="lower",
               vmin=0, vmax=vmax, extent=extent)
    ax2.axvline(bound_mb, color="blue", lw=1, linestyle="--", alpha=0.7)
    ax2.set_title(f"Contact: WT (left) | {best_label} (right)", fontweight="bold")
    ax2.set_xlabel("Genomic position (Mb)")
    ax2.set_ylabel("Genomic position (Mb)")

    fig.suptitle("Exp 3 — Synthetic Insulator Design at Unc5b", fontsize=13, fontweight="bold")
    os.makedirs("media", exist_ok=True)
    fig.tight_layout()
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {OUT_PNG}")
    print("Experiment 3 complete.")


if __name__ == "__main__":
    main()

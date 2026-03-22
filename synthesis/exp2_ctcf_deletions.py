#!/usr/bin/env python3
"""
Experiment 2 — CTCF deletion series.

Quantify insulation loss per CTCF site removed at Unc5b locus.

Variants tested:
  del_2ctcf  — delete 2 weakest CTCF peaks
  del_4ctcf  — delete all 4 weakest peaks
  flip_orient — flip strongest peak to divergent orientation

Output:
    media/synthesis_exp2_contact_panels.png
    media/synthesis_exp2_delta_insulation.png
"""

import os
import warnings
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as mgs
from dotenv import load_dotenv

# ── Constants ─────────────────────────────────────────────────────────────────
CHROM     = "chr10"
WIN_START = 60_240_000
WIN_END   = 61_288_576
CELL_TYPE = "EFO:0004038"
WT_CACHE  = "data/processed/synthesis_unc5b_wt_cache.npz"

CTCF_FWD  = "CCGCGAGGCGCAGG"
CTCF_REV  = "CCTGCGCCTCGCGG"

OUT_PANELS = "media/synthesis_exp2_contact_panels.png"
OUT_DELTA  = "media/synthesis_exp2_delta_insulation.png"


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
    """Load WT from cache. Returns None if missing or corrupt."""
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
        tf_meta   = wt_out.chip_tf.metadata
        ctcf_mask = tf_meta["transcription_factor"] == "CTCF"
        ctcf_signal = (
            wt_out.chip_tf.values[:, ctcf_mask].mean(axis=1)
            if ctcf_mask.sum() > 0
            else wt_out.chip_tf.values.mean(axis=1)
        )
    except Exception:
        ctcf_signal = wt_out.chip_tf.values.mean(axis=1)
    return dict(contact_mat=contact_mat, ctcf_signal=ctcf_signal,
                resolution=resolution, n_bins=n_bins)


def identify_peaks(ctcf_signal, n_peaks=4):
    """Return n_peaks peak indices sorted by signal ascending (weakest first)."""
    try:
        from scipy.signal import find_peaks
        threshold = np.nanpercentile(ctcf_signal, 60)
        peaks, _ = find_peaks(ctcf_signal, height=threshold, distance=3)
        if len(peaks) >= n_peaks:
            # Sort by signal ascending
            return peaks[np.argsort(ctcf_signal[peaks])][:n_peaks]
    except Exception:
        pass
    # Fallback
    sorted_idx = np.argsort(ctcf_signal)
    # Return top n_peaks (which are the strongest ones) but sorted by signal ascending
    return sorted_idx[-n_peaks:]


def genomic_pos(bin_idx, chip_tf_res):
    return WIN_START + int(bin_idx * chip_tf_res)


def validate_position(pos, label):
    if not (WIN_START <= pos < WIN_END):
        raise ValueError(
            f"Variant position {pos:,} is outside window "
            f"{WIN_START:,}–{WIN_END:,} for variant '{label}'"
        )


def predict_variant_safe(model, dna_client, interval, variant, label):
    """Predict variant; on failure log warning and return None."""
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

    # Load or recompute WT
    wt_data = load_wt_cache()
    if wt_data is None:
        wt_data = predict_wildtype(model, dna_client)

    contact_mat = wt_data["contact_mat"]
    ctcf_signal = wt_data["ctcf_signal"]
    resolution  = wt_data["resolution"]
    n_bins      = wt_data["n_bins"]

    chip_tf_res = (WIN_END - WIN_START) / len(ctcf_signal)
    peaks       = identify_peaks(ctcf_signal, n_peaks=4)   # weakest → strongest
    print(f"CTCF peaks (bin idx, weakest first): {peaks}")
    print(f"Chip-TF resolution: {chip_tf_res:,.0f} bp/bin")

    # Genomic positions
    p0       = genomic_pos(peaks[0], chip_tf_res)
    p1       = min(genomic_pos(peaks[1], chip_tf_res) + 500, WIN_END - 1)
    p3       = min(genomic_pos(peaks[3], chip_tf_res) + 500, WIN_END - 1)
    p_strong = genomic_pos(peaks[-1], chip_tf_res)

    # Only validate variant *start* positions (p1/p3 are span endpoints, not starts)
    for name, pos in [("del_2ctcf start", p0), ("flip_orient", p_strong)]:
        validate_position(pos, name)

    span_2 = p1 - p0
    span_4 = p3 - p0
    if span_2 < 1:
        span_2 = 500
    if span_4 < 1:
        span_4 = 1000

    interval = genome.Interval(chromosome=CHROM, start=WIN_START, end=WIN_END)

    variants = {
        "del_2ctcf": genome.Variant(
            chromosome=CHROM, position=p0 - 1,
            reference_bases="N" * span_2, alternate_bases="N",
            name="del_2ctcf",
        ),
        "del_4ctcf": genome.Variant(
            chromosome=CHROM, position=p0 - 1,
            reference_bases="N" * span_4, alternate_bases="N",
            name="del_4ctcf",
        ),
        "flip_orient": genome.Variant(
            chromosome=CHROM, position=p_strong - 1,
            reference_bases="N" * len(CTCF_REV), alternate_bases=CTCF_REV,
            name="flip_orient",
        ),
    }

    # Predict all variants
    mut_mats = {}
    for label, variant in variants.items():
        print(f"Predicting variant: {label} ...")
        mat = predict_variant_safe(model, dna_client, interval, variant, label)
        if mat is not None:
            mut_mats[label] = mat
        else:
            mut_mats[label] = contact_mat.copy()  # fallback to WT

    # Compute insulation
    ins_wt   = compute_insulation(contact_mat, resolution)
    bbin     = boundary_bin(resolution)
    ins_vals = {"WT": ins_wt}
    for lbl, mm in mut_mats.items():
        ins_vals[lbl] = compute_insulation(mm, resolution)

    # ── Figure 1: contact panels (2-row × 4-col) ──────────────────────────────
    labels_ordered = ["WT", "del_2ctcf", "del_4ctcf", "flip_orient"]
    titles         = ["WT", "Del 2 CTCF", "Del 4 CTCF", "Flip Orientation"]

    vmax_shared = np.nanpercentile(contact_mat[contact_mat > 0], 99)
    fig, axes   = plt.subplots(2, 4, figsize=(18, 9))

    mats_all = {"WT": contact_mat, **mut_mats}

    extent = [WIN_START / 1e6, WIN_END / 1e6, WIN_START / 1e6, WIN_END / 1e6]
    for col, (lbl, title) in enumerate(zip(labels_ordered, titles)):
        mat = mats_all[lbl]
        # Row 0: contact maps
        im0 = axes[0, col].imshow(mat, cmap="Reds", aspect="auto", origin="lower",
                                   vmin=0, vmax=vmax_shared, extent=extent)
        axes[0, col].set_title(title, fontsize=11)
        axes[0, col].tick_params(labelbottom=(col == 0), labelleft=(col == 0))

        # Row 1: difference maps (blank for WT)
        if col == 0:
            axes[1, col].axis("off")
        else:
            diff = mat - contact_mat
            vlim = np.nanpercentile(np.abs(diff), 98)
            axes[1, col].imshow(diff, cmap="RdBu_r", aspect="auto", origin="lower",
                                 vmin=-vlim, vmax=vlim, extent=extent)
            axes[1, col].set_title(f"Δ {title}", fontsize=11)
            axes[1, col].tick_params(labelbottom=True, labelleft=(col == 1))

    plt.colorbar(im0, ax=axes[0, :], fraction=0.01, pad=0.01, label="Contact freq.")
    fig.suptitle("Exp 2 — CTCF Deletion Series: Contact Maps", fontsize=13, fontweight="bold")
    os.makedirs("media", exist_ok=True)
    fig.savefig(OUT_PANELS, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Saved → {OUT_PANELS}")

    # ── Figure 2: Δ insulation bar chart ──────────────────────────────────────
    del_labels = ["del_2ctcf", "del_4ctcf", "flip_orient"]
    del_titles = ["Del 2 CTCF", "Del 4 CTCF", "Flip Orient"]
    delta_vals = [float(ins_vals[lbl][bbin] - ins_wt[bbin]) for lbl in del_labels]
    colors     = ["firebrick" if d < 0 else "steelblue" for d in delta_vals]

    fig2, ax = plt.subplots(figsize=(7, 4))
    bars = ax.barh(del_titles, delta_vals, color=colors)
    ax.axvline(0, color="black", lw=1)
    for bar, val in zip(bars, delta_vals):
        ax.text(val + (0.0001 if val >= 0 else -0.0001), bar.get_y() + bar.get_height() / 2,
                f"{val:+.4f}", va="center", ha="left" if val >= 0 else "right", fontsize=9)
    ax.set_xlabel("Δ Insulation score at boundary bin")
    ax.set_title("Exp 2 — Insulation Change at Unc5b Boundary", fontweight="bold")
    fig2.tight_layout()
    fig2.savefig(OUT_DELTA, dpi=150, bbox_inches="tight")
    plt.close(fig2)
    print(f"Saved → {OUT_DELTA}")
    print("Experiment 2 complete.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Experiment 1 — Baseline validation.

Confirm AlphaGenome predicts a known strong CTCF insulator cluster correctly
at the Unc5b locus (mm10 chr10:60,240,000–61,288,576).

Output:
    media/synthesis_exp1_baseline_unc5b.png
    data/processed/synthesis_unc5b_wt_cache.npz
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
WIN_END   = 61_288_576          # exactly 2^20 = 1,048,576 bp
CELL_TYPE = "EFO:0004038"       # Mouse ESC
WT_CACHE  = "data/processed/synthesis_unc5b_wt_cache.npz"
OUT_PNG   = "media/synthesis_exp1_baseline_unc5b.png"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_model():
    load_dotenv()
    from alphagenome.models import dna_client
    api_key = os.getenv("ALPHA_GENOME_API_KEY")
    if not api_key:
        raise EnvironmentError("ALPHA_GENOME_API_KEY not set in .env")
    model = dna_client.create(api_key)
    return model, dna_client


def compute_insulation(mat, resolution, window_bp=50_000):
    n   = mat.shape[0]
    win = max(2, int(window_bp / resolution))
    ins = np.full(n, np.nan)
    for i in range(win, n - win):
        ins[i] = mat[i - win:i, i:i + win].mean()
    return ins


def find_ctcf_peaks(ctcf_signal, n_peaks=6):
    """Return indices of CTCF peaks above 75th percentile."""
    try:
        from scipy.signal import find_peaks
        threshold = np.nanpercentile(ctcf_signal, 75)
        peaks, _ = find_peaks(ctcf_signal, height=threshold, distance=3)
        if len(peaks) >= 1:
            return peaks
    except Exception:
        pass
    # Fallback: top-N by signal
    return np.argsort(ctcf_signal)[-n_peaks:]


def predict_wildtype(model, dna_client):
    from alphagenome.data import genome
    interval = genome.Interval(chromosome=CHROM, start=WIN_START, end=WIN_END)
    organism = dna_client.Organism.MUS_MUSCULUS

    print(f"Predicting WT: {CHROM}:{WIN_START:,}–{WIN_END:,} ...")
    wt_out = model.predict_interval(
        interval=interval,
        organism=organism,
        requested_outputs={
            dna_client.OutputType.CONTACT_MAPS,
            dna_client.OutputType.CHIP_TF,
        },
        ontology_terms=[CELL_TYPE],
    )

    contact_mat = wt_out.contact_maps.values[:, :, 0]
    n_bins      = contact_mat.shape[0]
    resolution  = (WIN_END - WIN_START) / n_bins
    print(f"  Contact map: {n_bins}×{n_bins}, resolution ≈ {resolution:,.0f} bp/bin")

    # Extract CTCF signal
    try:
        tf_meta      = wt_out.chip_tf.metadata
        ctcf_mask    = tf_meta["transcription_factor"] == "CTCF"
        if ctcf_mask.sum() == 0:
            warnings.warn("No CTCF track found — using mean of all TF tracks.")
            ctcf_signal = wt_out.chip_tf.values.mean(axis=1)
        else:
            ctcf_signal = wt_out.chip_tf.values[:, ctcf_mask].mean(axis=1)
    except Exception as exc:
        warnings.warn(f"CHIP_TF parsing failed ({exc}); using all-TF mean.")
        ctcf_signal = wt_out.chip_tf.values.mean(axis=1)

    print(f"  CTCF signal shape: {ctcf_signal.shape}")
    return dict(
        contact_mat=contact_mat,
        ctcf_signal=ctcf_signal,
        resolution=resolution,
        n_bins=n_bins,
    )


def save_cache(wt_data):
    os.makedirs(os.path.dirname(WT_CACHE), exist_ok=True)
    np.savez_compressed(
        WT_CACHE,
        contact_mat=wt_data["contact_mat"],
        ctcf_signal=wt_data["ctcf_signal"],
        resolution=np.float64(wt_data["resolution"]),
        n_bins=np.int64(wt_data["n_bins"]),
    )
    print(f"WT cache saved → {WT_CACHE}")


def plot_baseline(wt_data, peak_indices):
    mat        = wt_data["contact_mat"]
    ctcf_sig   = wt_data["ctcf_signal"]
    resolution = wt_data["resolution"]
    n_bins     = wt_data["n_bins"]

    insulation = compute_insulation(mat, resolution)

    # X axis in Mb
    xbins_mb  = np.linspace(WIN_START / 1e6, WIN_END / 1e6, n_bins)
    xctcf_mb  = np.linspace(WIN_START / 1e6, WIN_END / 1e6, len(ctcf_sig))
    peak_mb   = xctcf_mb[peak_indices]

    vmax = np.nanpercentile(mat[mat > 0], 99)

    fig = plt.figure(figsize=(10, 12))
    gs  = mgs.GridSpec(3, 1, height_ratios=[2, 1, 1], hspace=0.08)

    # Panel 1 — contact heatmap
    ax0 = fig.add_subplot(gs[0])
    im  = ax0.imshow(
        mat, cmap="Reds", aspect="auto", origin="lower",
        vmin=0, vmax=vmax,
        extent=[WIN_START / 1e6, WIN_END / 1e6, WIN_START / 1e6, WIN_END / 1e6],
    )
    plt.colorbar(im, ax=ax0, fraction=0.03, pad=0.02, label="Contact freq.")
    for pm in peak_mb:
        ax0.axvline(pm, color="blue", lw=0.8, alpha=0.6)
        ax0.axhline(pm, color="blue", lw=0.8, alpha=0.6)
    ax0.set_ylabel("Genomic position (Mb)")
    ax0.set_title("Unc5b Locus — WT Contact Map (Mouse ESC)", fontsize=13, fontweight="bold")
    ax0.tick_params(labelbottom=False)

    # Panel 2 — CTCF signal
    ax1 = fig.add_subplot(gs[1], sharex=ax0)
    ax1.fill_between(xctcf_mb, ctcf_sig, alpha=0.6, color="steelblue", label="CTCF ChIP-seq")
    for pm in peak_mb:
        ax1.axvline(pm, color="navy", lw=1.0, linestyle="--", alpha=0.7)
    ax1.set_ylabel("CTCF signal")
    ax1.legend(loc="upper right", fontsize=8)
    ax1.tick_params(labelbottom=False)

    # Panel 3 — insulation score
    ax2 = fig.add_subplot(gs[2], sharex=ax0)
    ax2.plot(xbins_mb, insulation, color="darkgreen", lw=1.5, label="Insulation score")
    for pm in peak_mb:
        ax2.axvline(pm, color="darkgreen", lw=1.0, linestyle="--", alpha=0.7)
    ax2.set_ylabel("Insulation score")
    ax2.set_xlabel("Genomic position (Mb)")
    ax2.legend(loc="upper right", fontsize=8)

    os.makedirs("media", exist_ok=True)
    fig.savefig(OUT_PNG, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {OUT_PNG}")


# ── Main ──────────────────────────────────────────────────────────────────────

def main():
    model, dna_client = load_model()
    wt_data = predict_wildtype(model, dna_client)

    # Detect CTCF peaks
    peak_indices = find_ctcf_peaks(wt_data["ctcf_signal"], n_peaks=6)
    print(f"CTCF peaks detected at bins: {peak_indices}")

    save_cache(wt_data)
    plot_baseline(wt_data, peak_indices)
    print("Experiment 1 complete.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Stage 1 — Genomics context + protein sequence fetch.

Fetches the UNC5B (Q8BYU4) death domain sequence from UniProt and plots
the Unc5b locus CTCF signal from the existing WT cache (no new API call needed).

Outputs:
    data/processed/pipeline_target_sequence.fasta
    media/pipeline_stage1_expression.png
"""

import os
import json
import textwrap

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# ── Constants ──────────────────────────────────────────────────────────────────
UNC5B_ACCESSION    = "Q8K1S3"      # Mouse Netrin receptor UNC5B (canonical isoform)
DEATH_DOMAIN_START = 865          # 1-based residue index (UniProt annotation)
DEATH_DOMAIN_END   = 943
UNIPROT_BASE       = "https://rest.uniprot.org/uniprotkb"
WT_CACHE           = "data/processed/synthesis_unc5b_wt_cache.npz"
OUT_FASTA          = "data/processed/pipeline_target_sequence.fasta"
OUT_PNG            = "media/pipeline_stage1_expression.png"

# Genomic coordinates (mm10) matching WT cache
WIN_START     = 60_240_000
WIN_END       = 61_288_576
# UNC5B gene body approximate bounds (mm10 chr10)
GENE_START_MB = 60.50
GENE_END_MB   = 60.85


# ── Helpers ────────────────────────────────────────────────────────────────────

def fetch_unc5b_sequence(accession=UNC5B_ACCESSION,
                         start_res=DEATH_DOMAIN_START,
                         end_res=DEATH_DOMAIN_END) -> str:
    """Fetch protein sequence from UniProt; return death domain slice."""
    import urllib.request

    url = f"{UNIPROT_BASE}/{accession}.json"
    print(f"Fetching UniProt record: {url}")
    with urllib.request.urlopen(url, timeout=30) as resp:
        data = json.loads(resp.read().decode())

    full_seq = data["sequence"]["value"]
    print(f"  Full sequence length: {len(full_seq)} aa")

    # Try to auto-detect death domain from UniProt features
    detected_start, detected_end = start_res, end_res
    for feat in data.get("features", []):
        feat_type = feat.get("type", "")
        description = feat.get("description", "")
        if "death" in feat_type.lower() or "death" in description.lower():
            loc = feat.get("location", {})
            s = loc.get("start", {}).get("value")
            e = loc.get("end", {}).get("value")
            if s and e:
                detected_start, detected_end = int(s), int(e)
                print(f"  Death domain feature found: {detected_start}–{detected_end}")
                break
    else:
        print(f"  Death domain feature not found; using fallback {start_res}–{end_res}")

    domain_seq = full_seq[detected_start - 1 : detected_end]
    print(f"  Death domain sequence ({len(domain_seq)} aa): {domain_seq[:20]}...")
    return domain_seq


def save_fasta(sequence: str,
               accession: str = UNC5B_ACCESSION,
               start_res: int = DEATH_DOMAIN_START,
               end_res: int = DEATH_DOMAIN_END,
               out_path: str = OUT_FASTA) -> None:
    """Write sequence to FASTA with 60-char line wrap."""
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    header = f">{accession}|UNC5B_MOUSE death domain {start_res}-{end_res}"
    wrapped = "\n".join(textwrap.wrap(sequence, width=60))
    with open(out_path, "w") as f:
        f.write(header + "\n" + wrapped + "\n")
    print(f"FASTA saved → {out_path}")


def load_wt_cache(cache_path: str = WT_CACHE) -> dict | None:
    """Load WT cache. Returns None if missing or corrupt."""
    if not os.path.exists(cache_path):
        return None
    try:
        data = np.load(cache_path)
        return dict(
            ctcf_signal=data["ctcf_signal"],
            resolution=float(data["resolution"]),
            n_bins=int(data["n_bins"]),
        )
    except Exception as exc:
        print(f"Warning: WT cache could not be loaded ({exc}); skipping plot.")
        return None


def plot_locus_context(cache_data: dict | None, out_png: str = OUT_PNG) -> None:
    """Single-panel CTCF signal plot with UNC5B gene body span."""
    os.makedirs("media", exist_ok=True)

    fig, ax = plt.subplots(figsize=(10, 3.5))

    if cache_data is not None:
        ctcf_sig = cache_data["ctcf_signal"]
        x_mb = np.linspace(WIN_START / 1e6, WIN_END / 1e6, len(ctcf_sig))
        ax.fill_between(x_mb, ctcf_sig, alpha=0.65, color="steelblue", label="CTCF signal (AlphaGenome)")
        ax.plot(x_mb, ctcf_sig, color="steelblue", lw=0.8, alpha=0.8)
    else:
        ax.text(0.5, 0.5, "WT cache not found — run Stage 1 (baseline experiment) first.",
                transform=ax.transAxes, ha="center", va="center", color="#888",
                fontsize=11, fontstyle="italic")

    # UNC5B gene body span
    span = ax.axvspan(GENE_START_MB, GENE_END_MB, color="#FF6F00", alpha=0.18,
                      label=f"UNC5B gene body (~{GENE_START_MB}–{GENE_END_MB} Mb)")
    ax.axvline(GENE_START_MB, color="#FF6F00", lw=1.2, linestyle="--", alpha=0.7)
    ax.axvline(GENE_END_MB,   color="#FF6F00", lw=1.2, linestyle="--", alpha=0.7)

    # Annotation arrow
    ax.annotate(
        "UNC5B\n(death domain target)",
        xy=((GENE_START_MB + GENE_END_MB) / 2, ax.get_ylim()[1] if cache_data else 0.8),
        xytext=((GENE_START_MB + GENE_END_MB) / 2, ax.get_ylim()[1] * 1.15 if cache_data else 1.0),
        ha="center", va="bottom", fontsize=9, color="#BF360C",
        arrowprops=dict(arrowstyle="->", color="#BF360C", lw=1.2),
    )

    ax.set_xlabel("Genomic position (Mb)", fontsize=11)
    ax.set_ylabel("CTCF signal", fontsize=11)
    ax.set_title(
        f"Unc5b Locus — CTCF Signal Context (mm10 chr10, Mouse ESC)\n"
        f"AlphaGenome WT prediction: {WIN_START/1e6:.2f}–{WIN_END/1e6:.2f} Mb",
        fontsize=12, fontweight="bold",
    )
    ax.legend(loc="upper right", fontsize=9)
    ax.set_xlim(WIN_START / 1e6, WIN_END / 1e6)

    fig.tight_layout()
    fig.savefig(out_png, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved → {out_png}")


# ── Main ───────────────────────────────────────────────────────────────────────

def main():
    # Cache check
    if os.path.exists(OUT_FASTA):
        print(f"FASTA cache exists → {OUT_FASTA} (skip fetch; delete to rerun)")
        with open(OUT_FASTA) as f:
            lines = [l.strip() for l in f if not l.startswith(">")]
        seq = "".join(lines)
        print(f"  Cached sequence length: {len(seq)} aa")
    else:
        seq = fetch_unc5b_sequence()
        save_fasta(seq)

    cache = load_wt_cache()
    if cache is None:
        print(f"WT cache not found at {WT_CACHE}. Plot will show placeholder text.")
    plot_locus_context(cache)
    print("Stage 1 complete.")


if __name__ == "__main__":
    main()

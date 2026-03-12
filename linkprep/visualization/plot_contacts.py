"""
Visualize contact matrices, TADs, compartments, and loops.

Produces figures in the linkprep/figures/ directory.
"""

import os
import sys

import cooler
import matplotlib
matplotlib.use("Agg")
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
import config
from analysis.contact_analysis import (
    load_matrix,
    compute_insulation,
    call_tads,
    compute_compartments,
)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_region(region: str):
    chrom, rest = region.split(":")
    s, e = rest.split("-")
    return chrom, int(s.replace(",", "")), int(e.replace(",", ""))


def _region_label(region: str, clr: cooler.Cooler) -> str:
    chrom, start, end = _parse_region(region)
    return f"{chrom}:{start//1_000_000:.1f}–{end//1_000_000:.1f} Mb"


# ---------------------------------------------------------------------------
# Figure 1: Hi-C heatmap + insulation track
# ---------------------------------------------------------------------------

def plot_heatmap_with_insulation(clr: cooler.Cooler, region: str, out_path: str):
    chrom, start, end = _parse_region(region)
    mat = load_matrix(clr, region)

    ins_df = compute_insulation(clr, region,
                                windows_bp=config.INSULATION_WINDOWS)
    w_key  = f"w{config.INSULATION_WINDOWS[1] // 1000}k"
    ins    = ins_df[f"insulation_{w_key}"].values
    bounds = ins_df[f"is_boundary_{w_key}"].values.astype(bool)

    fig, axes = plt.subplots(
        2, 1, figsize=(9, 8),
        gridspec_kw={"height_ratios": [4, 1]},
    )

    # Heatmap
    vmax = np.nanpercentile(mat[mat > 0], 98) if (mat > 0).any() else 1.0
    axes[0].imshow(
        np.log1p(mat),
        cmap="Reds",
        origin="upper",
        aspect="auto",
        vmin=0,
        vmax=np.log1p(vmax),
    )
    axes[0].set_title(
        f"LinkPrep Contact Matrix — {_region_label(region, clr)}", fontsize=12
    )
    axes[0].set_xlabel("Genomic bin")
    axes[0].set_ylabel("Genomic bin")

    # Boundary lines
    for b_idx in np.where(bounds)[0]:
        axes[0].axhline(b_idx, color="cyan", lw=0.5, alpha=0.7)
        axes[0].axvline(b_idx, color="cyan", lw=0.5, alpha=0.7)

    # Insulation track
    axes[1].fill_between(range(len(ins)), ins, color="#1565C0", alpha=0.7)
    axes[1].axhline(0, color="black", lw=0.6, ls="--")
    for b_idx in np.where(bounds)[0]:
        axes[1].axvline(b_idx, color="red", lw=0.8, alpha=0.6)
    axes[1].set_ylabel("Insulation")
    axes[1].set_xlabel("Genomic bin")

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Figure 2: TAD structure overlay
# ---------------------------------------------------------------------------

def plot_tad_overlay(clr: cooler.Cooler, region: str, out_path: str):
    chrom, start, end = _parse_region(region)
    res = clr.binsize
    mat = load_matrix(clr, region)

    ins_df = compute_insulation(clr, region, windows_bp=[config.INSULATION_WINDOWS[1]])
    tad_df = call_tads(ins_df)

    fig, ax = plt.subplots(figsize=(8, 7))
    vmax = np.nanpercentile(mat[mat > 0], 98) if (mat > 0).any() else 1.0
    ax.imshow(np.log1p(mat), cmap="hot_r", origin="upper", aspect="auto",
              vmin=0, vmax=np.log1p(vmax))

    colors = plt.cm.Set2.colors
    for i, row in tad_df.iterrows():
        s_bin = (row["start"] - start) // res
        e_bin = (row["end"]   - start) // res
        color = colors[i % len(colors)]
        rect  = patches.Rectangle(
            (s_bin, s_bin), e_bin - s_bin, e_bin - s_bin,
            linewidth=1.5, edgecolor=color, facecolor="none",
        )
        ax.add_patch(rect)

    ax.set_title(f"TADs — {_region_label(region, clr)}\n"
                 f"({len(tad_df)} TADs called)", fontsize=11)
    ax.set_xlabel("Genomic bin")
    ax.set_ylabel("Genomic bin")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Figure 3: A/B compartments
# ---------------------------------------------------------------------------

def plot_compartments(clr: cooler.Cooler, region: str, out_path: str):
    chrom, start, end = _parse_region(region)
    mat   = load_matrix(clr, region)
    comp  = compute_compartments(clr, region)

    fig, axes = plt.subplots(
        2, 1, figsize=(10, 7),
        gridspec_kw={"height_ratios": [4, 1]},
    )

    vmax = np.nanpercentile(mat[mat > 0], 98) if (mat > 0).any() else 1.0
    axes[0].imshow(np.log1p(mat), cmap="Blues", origin="upper", aspect="auto",
                   vmin=0, vmax=np.log1p(vmax))
    axes[0].set_title(f"A/B Compartments — {_region_label(region, clr)}", fontsize=12)
    axes[0].set_xlabel("Genomic bin")
    axes[0].set_ylabel("Genomic bin")

    e1     = comp["E1"].values
    colors = np.where(e1 >= 0, "#E53935", "#1E88E5")
    axes[1].bar(range(len(e1)), e1, color=colors, width=1.0, linewidth=0)
    axes[1].axhline(0, color="black", lw=0.8)
    axes[1].set_ylabel("E1")
    axes[1].set_xlabel("Genomic bin")

    # Legend
    from matplotlib.patches import Patch
    axes[1].legend(
        handles=[Patch(color="#E53935", label="A compartment"),
                 Patch(color="#1E88E5", label="B compartment")],
        loc="upper right", fontsize=8,
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()
    print(f"  Saved → {out_path}")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def run_all(clr: cooler.Cooler):
    os.makedirs(config.FIGURES_DIR, exist_ok=True)
    for name, region in config.REGIONS.items():
        label = name.lower().replace(" ", "_")
        plot_heatmap_with_insulation(
            clr, region,
            os.path.join(config.FIGURES_DIR, f"contact_heatmap_{label}.png"),
        )
        plot_tad_overlay(
            clr, region,
            os.path.join(config.FIGURES_DIR, f"tads_{label}.png"),
        )
        plot_compartments(
            clr, region,
            os.path.join(config.FIGURES_DIR, f"compartments_{label}.png"),
        )


if __name__ == "__main__":
    if not os.path.exists(config.COOL_FILE):
        sys.exit(f"Cool file not found: {config.COOL_FILE}\n"
                 "Run: python sample_data/generate_sample_data.py")
    clr = cooler.Cooler(config.COOL_FILE)
    run_all(clr)

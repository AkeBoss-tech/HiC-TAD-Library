import cooler
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import os
import cooltools
from matplotlib.gridspec import GridSpec

from src.tad_boundaries import (
    parse_coordinates,
    make_view_df,
    compute_directionality_index,
    score_boundary_prominence,
    call_tad_intervals,
    multiscale_insulation,
    classify_boundary_persistence,
    boundary_pileup,
)

# --- Configuration (matches visualize_tads.py) ---
FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
RESOLUTION = 5000

regions = {
    "Sox11_Chr12": "chr12:26,000,000-28,000,000",
    "Mir9-2_Chr13": "chr13:83,500,000-84,500,000",
}

# Default insulation window for boundary detection at 5kb resolution
INSULATION_WINDOW = 50000


# ---------------------------------------------------------------------------
# 1. Directionality Index
# ---------------------------------------------------------------------------

def plot_directionality_index(file_path, region_name, coordinates, window_bp=500_000):
    print(f"Processing directionality index for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    di_df = compute_directionality_index(clr, coordinates, window_bp=window_bp)
    di = di_df['DI'].values
    x = np.arange(len(di))

    plt.figure(figsize=(12, 3))
    plt.fill_between(x, 0, di, where=(di > 0), color='red', alpha=0.5, interpolate=True)
    plt.fill_between(x, 0, di, where=(di < 0), color='blue', alpha=0.5, interpolate=True)
    plt.plot(x, di, color='black', lw=0.5)
    plt.axhline(0, color='grey', lw=0.5)
    plt.xlim(0, len(di))
    plt.ylabel("Directionality Index")
    plt.title(f"Directionality Index: {region_name}\n{coordinates} (window={window_bp // 1000}kb)")

    output_filename = f"media/{region_name}_directionality_index.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# 2. Boundary Strength
# ---------------------------------------------------------------------------

def plot_boundary_strength(file_path, region_name, coordinates, window_bp=None):
    if window_bp is None:
        window_bp = INSULATION_WINDOW
    print(f"Processing boundary strength for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    context = window_bp * 2
    fetch_start = max(0, start_bp - context)
    fetch_end = min(clr.chromsizes[chrom], end_bp + context)
    view_df = make_view_df(chrom, fetch_start, fetch_end)

    ins_table = cooltools.insulation(clr, [window_bp], ignore_diags=2, view_df=view_df)
    region_ins = ins_table[
        (ins_table['chrom'] == chrom) &
        (ins_table['start'] >= start_bp) &
        (ins_table['end'] <= end_bp)
    ].copy().reset_index(drop=True)

    scored = score_boundary_prominence(region_ins, window_bp)
    boundaries = scored[scored['boundary_class'].notna()].copy()
    boundaries = boundaries.sort_values('prominence', ascending=False).reset_index(drop=True)

    if len(boundaries) == 0:
        print(f"  No boundaries detected for {region_name}, skipping.")
        plt.close('all')
        return

    color_map = {'strong': '#d62728', 'weak': '#ff7f0e', 'sub_threshold': '#999999'}
    bar_colors = [color_map.get(c, '#999999') for c in boundaries['boundary_class']]

    plt.figure(figsize=(10, 5))
    x = np.arange(len(boundaries))
    plt.bar(x, boundaries['prominence'].values, color=bar_colors, edgecolor='none')
    plt.axhline(0.5, color='#d62728', ls='--', lw=1, label='Strong threshold')
    plt.axhline(0.2, color='#ff7f0e', ls='--', lw=1, label='Weak threshold')
    plt.xlabel("Boundary (ranked by prominence)")
    plt.ylabel("Prominence")
    plt.title(f"Boundary Strength: {region_name}\n{coordinates} (insulation window={window_bp // 1000}kb)")
    plt.legend(loc='upper right', fontsize='small')

    output_filename = f"media/{region_name}_boundary_strength.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# 3. Multi-scale Insulation Heatmap
# ---------------------------------------------------------------------------

def plot_multiscale_insulation(file_path, region_name, coordinates):
    print(f"Processing multi-scale insulation for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    ms_table, window_sizes = multiscale_insulation(clr, coordinates)

    # Build 2D matrix: rows = window sizes, cols = genomic bins
    n_bins = len(ms_table)
    n_windows = len(window_sizes)
    matrix = np.full((n_windows, n_bins), np.nan)

    for i, w in enumerate(window_sizes):
        col = f'log2_insulation_score_{w}'
        if col in ms_table.columns:
            matrix[i, :] = ms_table[col].values

    # Symmetric colormap centered at 0
    vmax = np.nanpercentile(np.abs(matrix[np.isfinite(matrix)]), 95) if np.any(np.isfinite(matrix)) else 1.0

    plt.figure(figsize=(12, 5))
    plt.imshow(
        matrix, aspect='auto', cmap='RdBu_r',
        vmin=-vmax, vmax=vmax,
        origin='lower', interpolation='none',
    )
    plt.colorbar(label='log2 Insulation Score')

    # Y-axis: window sizes in kb
    yticks = np.arange(n_windows)
    ylabels = [f'{w // 1000}kb' for w in window_sizes]
    plt.yticks(yticks, ylabels)
    plt.ylabel("Window Size")

    # X-axis: genomic position
    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    n_ticks = 5
    tick_positions = np.linspace(0, n_bins - 1, n_ticks).astype(int)
    tick_labels = [f'{(start_bp + t * RESOLUTION) / 1e6:.1f}Mb' for t in tick_positions]
    plt.xticks(tick_positions, tick_labels)
    plt.xlabel(f"Position on {chrom}")

    plt.title(f"Multi-Scale Insulation: {region_name}\n{coordinates}")

    output_filename = f"media/{region_name}_multiscale_insulation.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# 4. TAD Overlay on Contact Map
# ---------------------------------------------------------------------------

def plot_tad_overlay(file_path, region_name, coordinates, window_bp=None):
    if window_bp is None:
        window_bp = INSULATION_WINDOW
    print(f"Processing TAD overlay for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    chrom, start_bp, end_bp = parse_coordinates(coordinates)

    # Fetch Hi-C matrix
    matrix = clr.matrix(balance=True).fetch(coordinates)

    # Compute insulation and call TADs
    context = window_bp * 2
    fetch_start = max(0, start_bp - context)
    fetch_end = min(clr.chromsizes[chrom], end_bp + context)
    view_df = make_view_df(chrom, fetch_start, fetch_end)

    ins_table = cooltools.insulation(clr, [window_bp], ignore_diags=2, view_df=view_df)
    region_ins = ins_table[
        (ins_table['chrom'] == chrom) &
        (ins_table['start'] >= start_bp) &
        (ins_table['end'] <= end_bp)
    ].copy().reset_index(drop=True)

    scored = score_boundary_prominence(region_ins, window_bp)
    tads = call_tad_intervals(scored, clr, boundary_class_min='weak')
    # Filter TADs to target region
    tads = tads[
        (tads['chrom'] == chrom) &
        (tads['start'] >= start_bp) &
        (tads['end'] <= end_bp)
    ]

    # --- Plot ---
    fig = plt.figure(figsize=(10, 12))
    gs = GridSpec(2, 1, height_ratios=[3, 1], hspace=0.08)

    # Top: Heatmap with TAD triangles
    ax_heat = fig.add_subplot(gs[0])
    im = ax_heat.imshow(
        matrix, cmap='RdYlBu_r', interpolation='none',
        norm=colors.LogNorm(vmin=0.001, vmax=0.05),
    )

    # Draw TAD corner triangles
    for _, tad in tads.iterrows():
        s_bin = (tad['start'] - start_bp) / RESOLUTION
        e_bin = (tad['end'] - start_bp) / RESOLUTION
        # Draw L-shaped corner on the upper triangle
        ax_heat.plot([s_bin, s_bin], [s_bin, e_bin], color='black', lw=1.5)
        ax_heat.plot([s_bin, e_bin], [e_bin, e_bin], color='black', lw=1.5)

    ax_heat.set_title(f"TAD Overlay: {region_name}\n{coordinates}")
    fig.colorbar(im, ax=ax_heat, label="Contact Frequency", fraction=0.046, pad=0.04)

    # Bottom: Insulation track
    ax_ins = fig.add_subplot(gs[1])
    ins_col = f'log2_insulation_score_{window_bp}'
    x = np.arange(len(scored))
    ax_ins.plot(x, scored[ins_col].values, color='black', lw=1)
    ax_ins.axhline(0, color='grey', lw=0.5)
    ax_ins.set_ylabel(f"log2 Insulation\n({window_bp // 1000}kb)")
    ax_ins.set_xlim(0, len(scored))

    # Mark boundaries
    for cls, marker_color in [('strong', '#d62728'), ('weak', '#ff7f0e')]:
        bnds = scored[scored['boundary_class'] == cls]
        if len(bnds) > 0:
            ax_ins.scatter(
                bnds.index, bnds[ins_col].values,
                color=marker_color, s=30, zorder=5, label=cls,
            )
    ax_ins.legend(loc='lower left', fontsize='small')

    output_filename = f"media/{region_name}_tad_overlay.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# 5. Boundary Pileup
# ---------------------------------------------------------------------------

def plot_boundary_pileup(file_path, region_name, coordinates, window_bp=None, flank_bins=40):
    if window_bp is None:
        window_bp = INSULATION_WINDOW
    print(f"Processing boundary pileup for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    chrom, start_bp, end_bp = parse_coordinates(coordinates)

    # Compute insulation for the whole chromosome (pileup needs global boundaries)
    chrom_len = clr.chromsizes[chrom]
    view_df = make_view_df(chrom, 0, chrom_len)
    ins_table = cooltools.insulation(clr, [window_bp], ignore_diags=2, view_df=view_df)
    scored = score_boundary_prominence(ins_table, window_bp)

    pileup_matrix, n_boundaries = boundary_pileup(
        clr, scored, flank_bins=flank_bins, boundary_class_min='weak',
        normalize_by_expected=True,
    )

    if n_boundaries == 0:
        print(f"  No boundaries for pileup in {region_name}, skipping.")
        plt.close('all')
        return

    plt.figure(figsize=(6, 6))
    vmax = np.nanpercentile(pileup_matrix[np.isfinite(pileup_matrix)], 95) if np.any(np.isfinite(pileup_matrix)) else 2.0
    vmin = np.nanpercentile(pileup_matrix[np.isfinite(pileup_matrix)], 5) if np.any(np.isfinite(pileup_matrix)) else 0.5

    plt.imshow(
        pileup_matrix, cmap='RdYlBu_r', interpolation='none',
        vmin=vmin, vmax=vmax,
    )
    plt.colorbar(label='Obs / Expected')

    # Crosshairs at center
    center = flank_bins
    plt.axhline(center, color='black', ls='--', lw=0.8, alpha=0.5)
    plt.axvline(center, color='black', ls='--', lw=0.8, alpha=0.5)

    plt.title(f"Boundary Pileup (N={n_boundaries}): {region_name}\n{chrom} (insulation {window_bp // 1000}kb)")
    plt.xlabel("Bins from boundary")
    plt.ylabel("Bins from boundary")

    output_filename = f"media/{region_name}_boundary_pileup.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# 6. Combined Multi-Panel Figure
# ---------------------------------------------------------------------------

def plot_combined_boundary_panel(file_path, region_name, coordinates, window_bp=None):
    if window_bp is None:
        window_bp = INSULATION_WINDOW
    print(f"Processing combined boundary panel for {region_name}...")

    uri = f"{file_path}::resolutions/{RESOLUTION}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return

    chrom, start_bp, end_bp = parse_coordinates(coordinates)

    # --- Data ---
    matrix = clr.matrix(balance=True).fetch(coordinates)

    # Insulation + boundaries
    context = window_bp * 2
    fetch_start = max(0, start_bp - context)
    fetch_end = min(clr.chromsizes[chrom], end_bp + context)
    view_df = make_view_df(chrom, fetch_start, fetch_end)
    ins_table = cooltools.insulation(clr, [window_bp], ignore_diags=2, view_df=view_df)
    region_ins = ins_table[
        (ins_table['chrom'] == chrom) &
        (ins_table['start'] >= start_bp) &
        (ins_table['end'] <= end_bp)
    ].copy().reset_index(drop=True)
    scored = score_boundary_prominence(region_ins, window_bp)
    tads = call_tad_intervals(scored, clr, boundary_class_min='weak')
    tads = tads[
        (tads['chrom'] == chrom) &
        (tads['start'] >= start_bp) &
        (tads['end'] <= end_bp)
    ]

    # Directionality index
    di_df = compute_directionality_index(clr, coordinates, window_bp=window_bp)

    # --- Layout: 4 panels ---
    fig = plt.figure(figsize=(12, 14))
    gs = GridSpec(4, 1, height_ratios=[3, 1, 1, 0.3], hspace=0.08)
    n_bins = matrix.shape[0]

    # Panel 0: Heatmap + TAD triangles
    ax0 = fig.add_subplot(gs[0])
    im = ax0.imshow(
        matrix, cmap='RdYlBu_r', interpolation='none',
        norm=colors.LogNorm(vmin=0.001, vmax=0.05),
    )
    for _, tad in tads.iterrows():
        s = (tad['start'] - start_bp) / RESOLUTION
        e = (tad['end'] - start_bp) / RESOLUTION
        ax0.plot([s, s], [s, e], color='black', lw=1.5)
        ax0.plot([s, e], [e, e], color='black', lw=1.5)
    ax0.set_xlim(-0.5, n_bins - 0.5)
    ax0.set_ylim(n_bins - 0.5, -0.5)
    ax0.set_title(f"Combined Boundary Analysis: {region_name}\n{coordinates}")
    ax0.set_xticks([])
    fig.colorbar(im, ax=ax0, label="Contact Frequency", fraction=0.046, pad=0.04)

    # Panel 1: Insulation track
    ax1 = fig.add_subplot(gs[1])
    ins_col = f'log2_insulation_score_{window_bp}'
    x = np.arange(len(scored))
    ax1.plot(x, scored[ins_col].values, color='black', lw=1)
    ax1.axhline(0, color='grey', lw=0.5)
    ax1.set_ylabel(f"log2 Ins.\n({window_bp // 1000}kb)")
    ax1.set_xlim(0, n_bins)
    ax1.set_xticks([])
    for cls, mc in [('strong', '#d62728'), ('weak', '#ff7f0e')]:
        bnds = scored[scored['boundary_class'] == cls]
        if len(bnds) > 0:
            ax1.scatter(bnds.index, bnds[ins_col].values, color=mc, s=20, zorder=5, label=cls)
    ax1.legend(loc='lower left', fontsize='x-small')

    # Panel 2: DI track
    ax2 = fig.add_subplot(gs[2])
    di = di_df['DI'].values
    xdi = np.arange(len(di))
    ax2.fill_between(xdi, 0, di, where=(di > 0), color='red', alpha=0.4, interpolate=True)
    ax2.fill_between(xdi, 0, di, where=(di < 0), color='blue', alpha=0.4, interpolate=True)
    ax2.plot(xdi, di, color='black', lw=0.5)
    ax2.axhline(0, color='grey', lw=0.5)
    ax2.set_ylabel("DI")
    ax2.set_xlim(0, n_bins)
    ax2.set_xticks([])

    # Panel 3: Boundary annotation ticks
    ax3 = fig.add_subplot(gs[3])
    color_map = {'strong': '#d62728', 'weak': '#ff7f0e', 'sub_threshold': '#999999'}
    for idx, row in scored.iterrows():
        if row['boundary_class'] is not None:
            c = color_map.get(row['boundary_class'], '#999999')
            ax3.axvline(idx, color=c, lw=2, alpha=0.8)
    ax3.set_xlim(0, n_bins)
    ax3.set_yticks([])
    ax3.set_ylabel("Bounds", fontsize='small')

    # X-axis labels on bottom panel
    n_ticks = 6
    tick_positions = np.linspace(0, n_bins - 1, n_ticks).astype(int)
    tick_labels = [f'{(start_bp + t * RESOLUTION) / 1e6:.1f}Mb' for t in tick_positions]
    ax3.set_xticks(tick_positions)
    ax3.set_xticklabels(tick_labels, fontsize='small')
    ax3.set_xlabel(f"Position on {chrom}")

    output_filename = f"media/{region_name}_combined_boundary.png"
    plt.savefig(output_filename, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_filename}")
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if not os.path.exists(FILE_PATH):
        print(f"STOP: I can't find the file at {FILE_PATH}")
    else:
        for name, coords in regions.items():
            plot_directionality_index(FILE_PATH, name, coords)
            plot_boundary_strength(FILE_PATH, name, coords)
            plot_multiscale_insulation(FILE_PATH, name, coords)
            plot_tad_overlay(FILE_PATH, name, coords)
            plot_boundary_pileup(FILE_PATH, name, coords)
            plot_combined_boundary_panel(FILE_PATH, name, coords)

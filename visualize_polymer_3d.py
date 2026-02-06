"""
3-D polymer visualisation of Hi-C contact data.

Runs a lightweight bead-spring simulation (src/polymer_sim) on a genomic
region and renders the resulting structure coloured by A/B compartment,
insulation score, or TAD membership.

Two renderers are provided:
  * matplotlib  – static PNG, works everywhere
  * plotly      – interactive HTML, great in Jupyter

Usage (command-line)::

    python visualize_polymer_3d.py

Usage (library)::

    from visualize_polymer_3d import run_polymer_viz
    run_polymer_viz(file_path, "Sox11_Chr12", "chr12:26,000,000-28,000,000")
"""

import cooler
import cooltools
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from mpl_toolkits.mplot3d import Axes3D          # noqa: F401 (registers 3D projection)
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.colors as mcolors
import os

from src.tad_boundaries import (
    parse_coordinates,
    make_view_df,
    score_boundary_prominence,
    call_tad_intervals,
)
from src.polymer_sim import (
    contact_matrix_to_restraints,
    insulation_to_backbone_stiffness,
    simulate_polymer,
)


# --- Configuration (matches other visualize_*.py scripts) ---
FILE_PATH = "/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/data/raw/mouse_microc.mcool"
RESOLUTION = 5000
INSULATION_WINDOW = 50000

regions = {
    "Sox11_Chr12": "chr12:26,000,000-28,000,000",
    "Mir9-2_Chr13": "chr13:83,500,000-84,500,000",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _compute_eigenvector(clr, chrom: str) -> np.ndarray:
    """Return the E1 eigenvector for *chrom* (A/B compartment signal)."""
    view_df = pd.DataFrame({
        'chrom': [chrom], 'start': [0],
        'end': [clr.chromsizes[chrom]], 'name': [chrom],
    })
    _, eigvecs = cooltools.eigs_cis(clr, view_df=view_df, n_eigs=1)
    return eigvecs


def _region_insulation(clr, coordinates: str, window_bp: int) -> pd.DataFrame:
    """Insulation scores for *coordinates* with context padding."""
    chrom, start_bp, end_bp = parse_coordinates(coordinates)
    context = window_bp * 2
    fetch_start = max(0, start_bp - context)
    fetch_end = min(clr.chromsizes[chrom], end_bp + context)
    view_df = make_view_df(chrom, fetch_start, fetch_end)

    ins_table = cooltools.insulation(
        clr, [window_bp], ignore_diags=2, view_df=view_df,
    )
    region_ins = ins_table[
        (ins_table['chrom'] == chrom)
        & (ins_table['start'] >= start_bp)
        & (ins_table['end'] <= end_bp)
    ].copy().reset_index(drop=True)
    return region_ins


def _assign_tad_ids(scored: pd.DataFrame, clr, chrom: str,
                    start_bp: int, end_bp: int) -> np.ndarray:
    """Return per-bin TAD integer IDs (for coloring)."""
    tads = call_tad_intervals(scored, clr, boundary_class_min='weak')
    tads = tads[
        (tads['chrom'] == chrom)
        & (tads['start'] >= start_bp)
        & (tads['end'] <= end_bp)
    ]
    n_bins = len(scored)
    tad_ids = np.full(n_bins, -1, dtype=int)
    for _, tad in tads.iterrows():
        s = int((tad['start'] - start_bp) / clr.binsize)
        e = int((tad['end'] - start_bp) / clr.binsize)
        s, e = max(s, 0), min(e, n_bins)
        tad_ids[s:e] = tad['tad_id']
    return tad_ids


# ---------------------------------------------------------------------------
# Matplotlib renderer (static PNG)
# ---------------------------------------------------------------------------

def plot_polymer_matplotlib(
    coords: np.ndarray,
    bead_colors: np.ndarray,
    cmap: str = 'RdBu_r',
    title: str = '',
    output_path: str = 'media/polymer_3d.png',
    colorbar_label: str = 'E1 (A/B compartment)',
    elev: float = 25.0,
    azim: float = 135.0,
) -> None:
    """
    Render polymer as a 3-D scatter + backbone in matplotlib.

    Parameters
    ----------
    coords : np.ndarray
        (N, 3) bead positions.
    bead_colors : np.ndarray
        (N,) scalar values mapped through *cmap*.
    cmap : str
        Matplotlib colourmap name.
    title : str
        Figure title.
    output_path : str
        Where to save the PNG.
    colorbar_label : str
        Label for the colour bar.
    elev, azim : float
        Camera angles.
    """
    fig = plt.figure(figsize=(10, 10))
    ax = fig.add_subplot(111, projection='3d')

    # Backbone segments
    segments = [[coords[i], coords[i + 1]] for i in range(len(coords) - 1)]
    lc = Line3DCollection(segments, colors='#cccccc', linewidths=0.8, alpha=0.6)
    ax.add_collection3d(lc)

    # Beads
    vmax = np.nanpercentile(np.abs(bead_colors[np.isfinite(bead_colors)]), 95) \
        if np.any(np.isfinite(bead_colors)) else 1.0
    sc = ax.scatter(
        coords[:, 0], coords[:, 1], coords[:, 2],
        c=bead_colors, cmap=cmap, s=20, edgecolors='none',
        vmin=-vmax, vmax=vmax, depthshade=True,
    )
    fig.colorbar(sc, ax=ax, label=colorbar_label, shrink=0.6, pad=0.1)

    ax.set_title(title)
    ax.view_init(elev=elev, azim=azim)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()


def plot_polymer_with_heatmap(
    coords: np.ndarray,
    matrix: np.ndarray,
    bead_colors: np.ndarray,
    cmap_3d: str = 'RdBu_r',
    title: str = '',
    output_path: str = 'media/polymer_3d_panel.png',
    colorbar_label: str = 'E1 (A/B compartment)',
) -> None:
    """
    Side-by-side panel: 3-D polymer model (left) + contact heatmap (right).
    """
    fig = plt.figure(figsize=(16, 7))
    gs = GridSpec(1, 2, width_ratios=[1, 1], wspace=0.15)

    # --- Left: 3-D structure ---
    ax3d = fig.add_subplot(gs[0], projection='3d')
    segments = [[coords[i], coords[i + 1]] for i in range(len(coords) - 1)]
    lc = Line3DCollection(segments, colors='#cccccc', linewidths=0.8, alpha=0.6)
    ax3d.add_collection3d(lc)

    vmax = np.nanpercentile(np.abs(bead_colors[np.isfinite(bead_colors)]), 95) \
        if np.any(np.isfinite(bead_colors)) else 1.0
    sc = ax3d.scatter(
        coords[:, 0], coords[:, 1], coords[:, 2],
        c=bead_colors, cmap=cmap_3d, s=18, edgecolors='none',
        vmin=-vmax, vmax=vmax, depthshade=True,
    )
    fig.colorbar(sc, ax=ax3d, label=colorbar_label, shrink=0.5, pad=0.12)
    ax3d.set_title("3-D Polymer Model")
    ax3d.view_init(elev=25, azim=135)

    # --- Right: Contact heatmap ---
    ax_hm = fig.add_subplot(gs[1])
    im = ax_hm.imshow(
        matrix, cmap='RdYlBu_r', interpolation='none',
        norm=mcolors.LogNorm(vmin=0.001, vmax=0.05),
    )
    ax_hm.set_title("Hi-C Contact Map")
    fig.colorbar(im, ax=ax_hm, label='Contact Frequency', fraction=0.046, pad=0.04)

    fig.suptitle(title, fontsize=14, y=1.02)
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    print(f"Saved to {output_path}")
    plt.close()


# ---------------------------------------------------------------------------
# Plotly renderer (interactive HTML / Jupyter)
# ---------------------------------------------------------------------------

def plot_polymer_plotly(
    coords: np.ndarray,
    bead_colors: np.ndarray,
    cmap: str = 'RdBu',
    title: str = '',
    output_path: str = 'media/polymer_3d.html',
    colorbar_label: str = 'E1 (A/B compartment)',
):
    """
    Interactive 3-D polymer render using plotly.

    Returns the plotly Figure so it can be displayed inline in Jupyter via
    ``fig.show()``, or saved to *output_path* as standalone HTML.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        print("plotly not installed — skipping interactive render.")
        print("  Install with: pip install plotly")
        return None

    vmax = float(
        np.nanpercentile(np.abs(bead_colors[np.isfinite(bead_colors)]), 95)
    ) if np.any(np.isfinite(bead_colors)) else 1.0

    # Backbone trace
    backbone = go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        mode='lines',
        line=dict(color='grey', width=2),
        hoverinfo='skip',
        name='backbone',
    )

    # Bead trace
    beads = go.Scatter3d(
        x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
        mode='markers',
        marker=dict(
            size=3,
            color=bead_colors,
            colorscale=cmap,
            cmin=-vmax,
            cmax=vmax,
            colorbar=dict(title=colorbar_label),
            opacity=0.9,
        ),
        text=[f"bin {i}" for i in range(len(coords))],
        hoverinfo='text+x+y+z',
        name='beads',
    )

    fig = go.Figure(data=[backbone, beads])
    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title='X', yaxis_title='Y', zaxis_title='Z',
            aspectmode='data',
        ),
        width=800, height=700,
    )

    fig.write_html(output_path)
    print(f"Saved interactive HTML to {output_path}")
    return fig


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------

def run_polymer_viz(
    file_path: str,
    region_name: str,
    coordinates: str,
    resolution: int = RESOLUTION,
    insulation_window: int = INSULATION_WINDOW,
    n_steps: int = 5000,
    dt: float = 0.005,
    friction: float = 1.0,
    temperature: float = 1.0,
    contact_threshold_quantile: float = 0.70,
    seed: int = 42,
    render_plotly: bool = True,
) -> np.ndarray:
    """
    Full pipeline: load data → simulate → render.

    Parameters
    ----------
    file_path : str
        Path to .mcool file.
    region_name : str
        Human-readable label for the region.
    coordinates : str
        Genomic coordinates, e.g. "chr12:26,000,000-28,000,000".
    resolution : int
        Bin size in bp.
    insulation_window : int
        Insulation window for backbone stiffness modulation.
    n_steps, dt, friction, temperature : float
        Simulation parameters forwarded to ``simulate_polymer``.
    contact_threshold_quantile : float
        Quantile filter for Hi-C restraints.
    seed : int
        Random seed.
    render_plotly : bool
        If True, also produce an interactive HTML file via plotly.

    Returns
    -------
    np.ndarray
        Final bead coordinates (N, 3).
    """
    print(f"Processing 3-D polymer for {region_name} ({coordinates})...")

    uri = f"{file_path}::resolutions/{resolution}"
    try:
        clr = cooler.Cooler(uri)
    except Exception as e:
        print(f"Error: Could not load file. {e}")
        return np.array([])

    chrom, start_bp, end_bp = parse_coordinates(coordinates)

    # --- 1. Contact matrix ---
    print("  Fetching contact matrix...")
    matrix = clr.matrix(balance=True).fetch(coordinates)
    n_beads = matrix.shape[0]
    print(f"  Region has {n_beads} bins at {resolution}bp resolution.")

    # --- 2. Insulation → backbone stiffness ---
    print("  Computing insulation scores...")
    region_ins = _region_insulation(clr, coordinates, insulation_window)
    ins_col = f'log2_insulation_score_{insulation_window}'
    ins_values = region_ins[ins_col].values[:n_beads] if ins_col in region_ins.columns else None
    backbone_k = insulation_to_backbone_stiffness(ins_values) if ins_values is not None else None

    # --- 3. A/B compartment eigenvector ---
    print("  Computing A/B compartment eigenvector...")
    # Use 100kb resolution for compartment calls (standard practice)
    comp_res = 100_000
    comp_uri = f"{file_path}::resolutions/{comp_res}"
    try:
        clr_comp = cooler.Cooler(comp_uri)
        eigvecs = _compute_eigenvector(clr_comp, chrom)
        e1_region = eigvecs[
            (eigvecs['chrom'] == chrom)
            & (eigvecs['start'] >= start_bp)
            & (eigvecs['end'] <= end_bp)
        ]['E1'].values
        # Upsample to match n_beads (100kb → 5kb = repeat each value 20x)
        scale = resolution
        comp_scale = comp_res
        if len(e1_region) > 0:
            e1_upsampled = np.repeat(e1_region, comp_scale // scale)[:n_beads]
            # Pad if rounding mismatch
            if len(e1_upsampled) < n_beads:
                e1_upsampled = np.pad(e1_upsampled, (0, n_beads - len(e1_upsampled)),
                                      constant_values=0.0)
        else:
            e1_upsampled = np.zeros(n_beads)
    except Exception as e:
        print(f"  Warning: could not compute eigenvector ({e}), using insulation for colour.")
        e1_upsampled = ins_values if ins_values is not None else np.zeros(n_beads)

    # --- 4. Simulate ---
    print(f"  Building restraints (threshold quantile={contact_threshold_quantile})...")
    restraints = contact_matrix_to_restraints(
        matrix, contact_threshold_quantile=contact_threshold_quantile,
    )
    print(f"  {len(restraints)} Hi-C restraints + {n_beads - 1} backbone springs.")

    print(f"  Running simulation ({n_steps} steps, dt={dt}, friction={friction})...")
    coords = simulate_polymer(
        n_beads, restraints, backbone_k=backbone_k,
        n_steps=n_steps, dt=dt, friction=friction,
        temperature=temperature, seed=seed,
    )
    print("  Simulation complete.")

    # --- 5. Render ---
    # Standalone 3-D
    plot_polymer_matplotlib(
        coords, e1_upsampled,
        title=f"3-D Polymer: {region_name}\n{coordinates}",
        output_path=f"media/{region_name}_polymer_3d.png",
    )

    # Side-by-side panel
    plot_polymer_with_heatmap(
        coords, matrix, e1_upsampled,
        title=f"{region_name} — {coordinates}",
        output_path=f"media/{region_name}_polymer_3d_panel.png",
    )

    # Interactive plotly
    if render_plotly:
        plot_polymer_plotly(
            coords, e1_upsampled,
            title=f"3-D Polymer: {region_name} ({coordinates})",
            output_path=f"media/{region_name}_polymer_3d.html",
        )

    return coords


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    if not os.path.exists(FILE_PATH):
        print(f"STOP: I can't find the file at {FILE_PATH}")
    else:
        for name, coords in regions.items():
            run_polymer_viz(FILE_PATH, name, coords)

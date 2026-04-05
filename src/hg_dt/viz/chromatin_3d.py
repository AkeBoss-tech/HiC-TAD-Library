"""
3D chromatin conformation visualization.

Runs the polymer simulation on ref and mut contact maps, interpolates
between the two equilibrium structures, and renders an animated Plotly
scatter3d showing how the chromatin reorganizes after the edit.
"""

import numpy as np
import streamlit as st
import plotly.graph_objects as go
from typing import Optional


def _downsample_matrix(matrix: np.ndarray, target: int = 30) -> np.ndarray:
    """Block-average a contact matrix down to at most `target` bins."""
    n = matrix.shape[0]
    if n <= target:
        return matrix
    block = n // target
    n_out = n // block
    m = np.zeros((n_out, n_out))
    for i in range(n_out):
        for j in range(n_out):
            m[i, j] = np.nanmean(matrix[i*block:(i+1)*block, j*block:(j+1)*block])
    return m


@st.cache_data(show_spinner="Simulating 3D chromatin conformation…")
def run_chromatin_sim(
    ref_hic: np.ndarray,
    mut_hic: np.ndarray,
    n_frames: int = 40,
    n_steps: int = 2000,
) -> tuple:
    """
    Run polymer simulation on ref and mut contact maps, return interpolated
    trajectory frames.

    Returns
    -------
    frames : np.ndarray  shape (n_frames+1, n_beads, 3)
    ref_coords : np.ndarray  shape (n_beads, 3)
    mut_coords : np.ndarray  shape (n_beads, 3)
    """
    from src.polymer_sim import contact_matrix_to_restraints, simulate_polymer

    ref_m = _downsample_matrix(ref_hic, target=30)
    mut_m = _downsample_matrix(mut_hic, target=30)
    n = ref_m.shape[0]

    ref_restraints = contact_matrix_to_restraints(ref_m)
    mut_restraints = contact_matrix_to_restraints(mut_m)

    ref_coords = simulate_polymer(n, ref_restraints, n_steps=n_steps, seed=42)
    mut_coords = simulate_polymer(n, mut_restraints, n_steps=n_steps, seed=42)

    # Align mut_coords to ref_coords via centroid + Kabsch rotation
    ref_c = ref_coords - ref_coords.mean(axis=0)
    mut_c = mut_coords - mut_coords.mean(axis=0)
    H = mut_c.T @ ref_c
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    mut_aligned = (mut_c @ R.T) + ref_coords.mean(axis=0)

    frames = np.array([
        (1 - t) * ref_coords + t * mut_aligned
        for t in np.linspace(0, 1, n_frames + 1)
    ])
    return frames, ref_coords, mut_aligned


def render_chromatin_animation(
    ref_hic: np.ndarray,
    mut_hic: np.ndarray,
    edit_start: int,
    edit_end: int,
    genomic_start: int,
    genomic_end: int,
    n_frames: int = 40,
):
    """
    Render an animated 3D polymer showing the transition from reference
    to mutant chromatin conformation.

    Beads are colored by genomic position (rainbow). Beads that overlap
    the edit region are highlighted in red.
    """
    frames, ref_coords, mut_coords = run_chromatin_sim(ref_hic, mut_hic, n_frames)
    n_beads = ref_coords.shape[0]

    # Genomic position of each bead
    bead_positions = np.linspace(genomic_start, genomic_end, n_beads)
    in_edit = (bead_positions >= edit_start) & (bead_positions <= edit_end)

    # Color: rainbow for all beads, bright red for edit region
    colors = np.linspace(0, 1, n_beads)

    def _make_trace(coords, frame_idx):
        # Lines (backbone)
        line_x, line_y, line_z = [], [], []
        for i in range(n_beads - 1):
            line_x += [coords[i, 0], coords[i+1, 0], None]
            line_y += [coords[i, 1], coords[i+1, 1], None]
            line_z += [coords[i, 2], coords[i+1, 2], None]

        marker_colors = [
            'rgba(255,50,50,1)' if in_edit[i] else
            f'hsl({int(300 * (1 - colors[i]))}, 80%, 55%)'
            for i in range(n_beads)
        ]
        marker_sizes = [10 if in_edit[i] else 6 for i in range(n_beads)]

        return [
            go.Scatter3d(
                x=line_x, y=line_y, z=line_z,
                mode='lines',
                line=dict(color='rgba(200,200,200,0.4)', width=3),
                hoverinfo='none',
                showlegend=False,
            ),
            go.Scatter3d(
                x=coords[:, 0], y=coords[:, 1], z=coords[:, 2],
                mode='markers',
                marker=dict(
                    size=marker_sizes,
                    color=marker_colors,
                    opacity=0.92,
                    line=dict(width=0.5, color='rgba(0,0,0,0.3)'),
                ),
                text=[
                    f"Bead {i}<br>~{int(bead_positions[i]):,} bp"
                    + (" ← edit region" if in_edit[i] else "")
                    for i in range(n_beads)
                ],
                hoverinfo='text',
                showlegend=False,
            ),
        ]

    # Build animation frames
    plotly_frames = [
        go.Frame(data=_make_trace(frames[i], i), name=str(i))
        for i in range(len(frames))
    ]

    # Axis ranges
    all_pts = frames.reshape(-1, 3)
    pad = (all_pts.max() - all_pts.min()) * 0.1 + 0.1
    rng = [float(all_pts.min() - pad), float(all_pts.max() + pad)]
    axis = dict(range=rng, showgrid=False, zeroline=False,
                showticklabels=False, showspikes=False, backgroundcolor='#0d1117')

    fig = go.Figure(
        data=_make_trace(frames[0], 0),
        frames=plotly_frames,
        layout=go.Layout(
            scene=dict(
                xaxis=axis, yaxis=axis, zaxis=axis,
                bgcolor='#0d1117',
                camera=dict(eye=dict(x=1.4, y=1.4, z=0.8)),
            ),
            paper_bgcolor='#0d1117',
            plot_bgcolor='#0d1117',
            font=dict(color='#e6edf3'),
            height=520,
            margin=dict(l=0, r=0, t=40, b=0),
            title=dict(
                text='3D Chromatin Conformation: Reference → Mutant',
                font=dict(size=14, color='#e6edf3'),
            ),
            updatemenus=[dict(
                type='buttons',
                showactive=False,
                y=0.02, x=0.5, xanchor='center', yanchor='bottom',
                buttons=[
                    dict(
                        label='▶  Play transition',
                        method='animate',
                        args=[None, dict(
                            frame=dict(duration=80, redraw=True),
                            fromcurrent=True,
                            transition=dict(duration=40),
                        )],
                    ),
                    dict(
                        label='⏸  Pause',
                        method='animate',
                        args=[[None], dict(
                            frame=dict(duration=0, redraw=False),
                            mode='immediate',
                            transition=dict(duration=0),
                        )],
                    ),
                ],
                bgcolor='#21262d',
                font=dict(color='#e6edf3'),
                bordercolor='#30363d',
            )],
            sliders=[dict(
                active=0,
                currentvalue=dict(
                    prefix='Frame: ',
                    visible=True,
                    xanchor='right',
                    font=dict(color='#8b949e'),
                ),
                pad=dict(t=50, b=10),
                bgcolor='#21262d',
                bordercolor='#30363d',
                tickcolor='#8b949e',
                font=dict(color='#8b949e'),
                steps=[
                    dict(
                        args=[[f.name], dict(
                            frame=dict(duration=80, redraw=True),
                            mode='immediate',
                            transition=dict(duration=0),
                        )],
                        method='animate',
                        label='',
                    )
                    for f in plotly_frames
                ],
            )],
        ),
    )

    st.plotly_chart(fig, use_container_width=True)

    # Compact legend
    lcol1, lcol2, lcol3 = st.columns(3)
    with lcol1:
        st.caption("🔵→🔴 Bead color = genomic position (5′ → 3′)")
    with lcol2:
        st.caption("🔴 Bright red beads = edit region")
    with lcol3:
        st.caption("Gray line = chromatin backbone")

import streamlit as st
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
from PIL import Image

# Import the new modules
from src.hg_dt.viz.tracks_plotter import plot_tracks
from src.hg_dt.viz.hic_plotter import plot_hic_triangle
from src.hg_dt.viz.protein_viz import render_protein_comparison
from src.hg_dt.analyze.attribution import generate_mechanistic_insight

st.set_page_config(page_title="HG-DT Causal Dashboard", layout="wide")

st.title("HG-DT: Visual Interpretation & Causal Dashboard")

# Mock data generation for demonstration
@st.cache_data
def generate_mock_data(mod_type: str, gene: str):
    # Tracks: 100 bins
    ref_track = np.random.normal(5, 1, 100)
    mut_track = ref_track.copy()
    if mod_type == "Deletion":
        # Simulate loss in middle bins
        mut_track[40:60] = mut_track[40:60] * 0.2

    # Contact Maps: 50x50
    ref_hic = np.random.rand(50, 50)
    ref_hic = (ref_hic + ref_hic.T) / 2
    np.fill_diagonal(ref_hic, 1.0)

    mut_hic = ref_hic.copy()
    if mod_type == "Deletion":
        # Disrupted loop
        mut_hic[10:20, 30:40] = mut_hic[10:20, 30:40] * 0.3
        mut_hic[30:40, 10:20] = mut_hic[10:20, 30:40]

    mock_pdb = ""

    return ref_track, mut_track, ref_hic, mut_hic, mock_pdb, mock_pdb

@st.cache_data
def generate_accessibility_over_time(gene: str, mod_type: str, time_steps: int = 50):
    """Simulate gene accessibility changing over time during structural transition."""
    t = np.linspace(0, 10, time_steps)
    base_accessibility = 100.0

    if mod_type == "Deletion":
        # Sigmoidal drop in accessibility
        accessibility = base_accessibility - 60 * (1 / (1 + np.exp(-t + 5))) + np.random.normal(0, 2, time_steps)
    elif mod_type == "Insertion":
        # Sigmoidal increase
        accessibility = base_accessibility + 80 * (1 / (1 + np.exp(-t + 5))) + np.random.normal(0, 2, time_steps)
    else:
        accessibility = np.full(time_steps, base_accessibility) + np.random.normal(0, 2, time_steps)

    return t, accessibility

# Tabs definition
tab1, tab2, tab3, tab4, tab5, tab6 = st.tabs([
    "Specification",
    "Genome Tracks",
    "3D Organization",
    "Protein Structure",
    "Trajectory Animation",
    "Mechanistic Attribution"
])

if "mod_type" not in st.session_state:
    st.session_state.mod_type = "None"
if "gene" not in st.session_state:
    st.session_state.gene = "TAL1"
if "analyzed" not in st.session_state:
    st.session_state.analyzed = False

with tab1:
    st.header("Input DNA Modification")

    col1, col2 = st.columns(2)
    with col1:
        gene = st.selectbox(
            "Select Target Gene:",
            ["TAL1", "OCT4", "NANOG", "SOX2", "Mef2c"]
        )
    with col2:
        mod_type = st.selectbox(
            "Select Modification Type:",
            ["None", "Deletion", "Insertion", "Duplication"]
        )

    chrom = st.text_input("Chromosome", "chr1")
    locus = st.text_input("Locus", "47200000-47250000")

    if st.button("Run Analysis"):
        st.session_state.gene = gene
        st.session_state.mod_type = mod_type
        st.session_state.analyzed = True
        st.success(f"Analysis triggered for {gene} ({mod_type}) at {chrom}:{locus}")

if st.session_state.analyzed and st.session_state.mod_type != "None":
    ref_track, mut_track, ref_hic, mut_hic, ref_pdb, mut_pdb = generate_mock_data(
        st.session_state.mod_type, st.session_state.gene
    )

    with tab2:
        st.header("1D Genome Tracks (Ref vs. Mut)")
        st.write(f"Linear browser views for accessibility and expression deltas for **{st.session_state.gene}**.")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            track_img_path = plot_tracks(ref_track, mut_track, tmp.name, title=f"1D Tracks: {st.session_state.gene} {st.session_state.mod_type}")
        st.image(Image.open(track_img_path), use_column_width=True)

    with tab3:
        st.header("3D Organization (Hi-C Heatmaps)")
        st.write("Triangular/Matrix contact map deltas + loop reorganization.")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            hic_img_path = plot_hic_triangle(ref_hic, mut_hic, tmp.name, title=f"3D Contact Map: {st.session_state.gene} {st.session_state.mod_type}")
        st.image(Image.open(hic_img_path), use_column_width=True)

    with tab4:
        st.header("Protein Structure")
        st.write("View 3D folding (Ref vs. Mut).")
        render_protein_comparison(ref_pdb, mut_pdb)

    with tab5:
        st.header("Trajectory Animation")
        st.write(f"Molecular simulation trajectory: tracking **{st.session_state.gene}** accessibility over time as structure transitions from reference to mutant state.")

        t, accessibility = generate_accessibility_over_time(st.session_state.gene, st.session_state.mod_type)

        fig, ax = plt.subplots(figsize=(10, 4))
        ax.plot(t, accessibility, marker='o', linestyle='-', color='purple')
        ax.set_title(f"{st.session_state.gene} Accessibility Over Simulation Time")
        ax.set_xlabel("Simulation Step (Time)")
        ax.set_ylabel("Accessibility Signal")
        ax.grid(True, alpha=0.3)

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            fig.tight_layout()
            fig.savefig(tmp.name, dpi=150)
            plt.close(fig)
            st.image(Image.open(tmp.name), use_column_width=True)

        st.info("The plot above simulates the real-time changes in accessibility as the chromatin refolds post-modification.")

    with tab6:
        st.header("Mechanistic Attribution")

        mod_details = {"type": st.session_state.mod_type.lower(), "target": f"enhancer element near {st.session_state.gene}"}

        if st.session_state.mod_type == "Deletion":
            delta_stats = {
                "loop_weakened": True,
                "accessibility_drop": 0.28 if st.session_state.gene == "TAL1" else 0.45,
                "expression_drop": 0.35 if st.session_state.gene == "TAL1" else 0.50
            }
        elif st.session_state.mod_type == "Insertion":
            delta_stats = {
                "loop_strengthened": True,
                "accessibility_drop": -0.50, # Negative drop = gain
                "expression_drop": -0.80
            }
        else:
            delta_stats = {
                "loop_weakened": False,
                "accessibility_drop": 0.0,
                "expression_drop": 0.0
            }

        insight = generate_mechanistic_insight(mod_details, delta_stats)

        st.info("### Mechanistic Summary")
        st.success(insight)

elif st.session_state.analyzed:
    st.warning("Please select a valid modification from the Specification tab.")

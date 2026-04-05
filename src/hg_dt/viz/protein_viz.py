import py3Dmol
import streamlit as st
from stmol import showmol

def show_protein_3d(pdb_data: str, width: int = 800, height: int = 400):
    """
    Render a 3D protein structure using py3Dmol in Streamlit.

    Args:
        pdb_data: String containing the PDB file data.
        width: Width of the viewer.
        height: Height of the viewer.
    """
    view = py3Dmol.view(width=width, height=height)
    view.addModel(pdb_data, "pdb")
    view.setStyle({'cartoon': {'color': 'spectrum'}})
    view.zoomTo()
    showmol(view, height=height, width=width)

def render_protein_comparison(ref_pdb: str, mut_pdb: str):
    """
    Render a side-by-side comparison of Ref and Mut proteins in Streamlit.
    """
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Reference Structure")
        if ref_pdb:
            show_protein_3d(ref_pdb, width=350, height=350)
        else:
            st.info("Reference PDB not provided.")

    with col2:
        st.subheader("Mutant Structure")
        if mut_pdb:
            show_protein_3d(mut_pdb, width=350, height=350)
        else:
            st.info("Mutant PDB not provided.")

import py3Dmol
import streamlit as st
from stmol import showmol

def show_protein_3d(pdb_data: str, width: int = 400, height: int = 400, animate=False):
    """
    Render a 3D protein structure using py3Dmol in Streamlit.
    """
    view = py3Dmol.view(width=width, height=height)
    
    if "MODEL" in pdb_data:
        view.addModelsAsFrames(pdb_data, "pdb")
    else:
        view.addModel(pdb_data, "pdb")
        
    # Styles for simple mock coordinates (CA atoms)
    view.setStyle({'stick': {}, 'sphere': {'radius': 1.0}})
    view.setBackgroundColor('white')
    view.zoomTo()
    
    if animate:
        view.animate({'loop': 'forward', 'step': 2})
        view.spin(True)
        
    showmol(view, height=height, width=width)

def render_protein_comparison(ref_pdb: str, mut_pdb: str, animate=True):
    """
    Render a side-by-side comparison of Ref and Mut proteins in Streamlit.
    """
    col1, col2 = st.columns(2)

    with col1:
        st.subheader("Reference Structure (Dynamic)")
        if ref_pdb:
            show_protein_3d(ref_pdb, width=400, height=400, animate=animate)
        else:
            st.info("Reference PDB not provided.")

    with col2:
        st.subheader("Mutant Structure (Dynamic)")
        if mut_pdb:
            show_protein_3d(mut_pdb, width=400, height=400, animate=animate)
        else:
            st.info("Mutant PDB not provided.")

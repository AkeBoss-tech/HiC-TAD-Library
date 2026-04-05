import py3Dmol
import streamlit as st
from stmol import showmol

# Colors
_REF_COLOR = "#4a9eff"   # blue  — reference
_MUT_COLOR = "#ff6b35"   # orange — mutant
_BG_COLOR  = "#0d1117"   # dark background


def _has_full_backbone(pdb_text: str) -> bool:
    """Return True if the PDB has N/C/O atoms (real ESMFold) vs CA-only mock."""
    for line in pdb_text.splitlines():
        if line.startswith("ATOM") and " N " in line:
            return True
    return False


def show_protein_3d(pdb_data: str, width: int = 400, height: int = 400, animate=False,
                    color: str = _REF_COLOR):
    """
    Render a single protein structure using py3Dmol in Streamlit.
    Uses cartoon for real ESMFold output, spheres+sticks for CA-only mock.
    """
    view = py3Dmol.view(width=width, height=height)

    if "MODEL" in pdb_data:
        view.addModelsAsFrames(pdb_data, "pdb")
    else:
        view.addModel(pdb_data, "pdb")

    if _has_full_backbone(pdb_data):
        view.setStyle({'cartoon': {'color': 'spectrum', 'style': 'oval'}})
    else:
        view.setStyle({'sphere': {'radius': 0.5, 'color': color},
                       'stick':  {'radius': 0.15, 'color': color}})

    view.setBackgroundColor(_BG_COLOR)
    view.zoomTo()

    if animate:
        view.animate({'loop': 'forward', 'step': 1, 'interval': 80})
        view.spin(True)

    showmol(view, height=height, width=width)


def render_protein_overlay(ref_pdb: str, mut_pdb: str, width: int = 660, height: int = 440):
    """
    Render reference (blue) and mutant (orange) proteins overlaid in one viewer.

    Both structures are loaded into the same py3Dmol scene and colored differently
    so structural changes are immediately visible. Real ESMFold output gets
    cartoon ribbons; the CA-only mock gets colored spheres.
    """
    if not ref_pdb and not mut_pdb:
        st.info("No protein structures available — run with a real sequence to generate them.")
        return

    view = py3Dmol.view(width=width, height=height)
    model_idx = 0

    if ref_pdb:
        view.addModel(ref_pdb, 'pdb')
        if _has_full_backbone(ref_pdb):
            view.setStyle(
                {'model': model_idx},
                {'cartoon': {'color': _REF_COLOR, 'opacity': 0.85}},
            )
        else:
            view.setStyle(
                {'model': model_idx},
                {'sphere': {'radius': 0.45, 'color': _REF_COLOR},
                 'stick':  {'radius': 0.12, 'color': _REF_COLOR}},
            )
        model_idx += 1

    if mut_pdb:
        view.addModel(mut_pdb, 'pdb')
        if _has_full_backbone(mut_pdb):
            view.setStyle(
                {'model': model_idx},
                {'cartoon': {'color': _MUT_COLOR, 'opacity': 0.85}},
            )
        else:
            view.setStyle(
                {'model': model_idx},
                {'sphere': {'radius': 0.45, 'color': _MUT_COLOR},
                 'stick':  {'radius': 0.12, 'color': _MUT_COLOR}},
            )

    view.setBackgroundColor(_BG_COLOR)
    view.zoomTo()
    view.spin(True)

    showmol(view, height=height, width=width)

    # Inline legend
    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"<span style='color:{_REF_COLOR}; font-weight:600;'>■ Reference</span>",
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            f"<span style='color:{_MUT_COLOR}; font-weight:600;'>■ Mutant</span>",
            unsafe_allow_html=True,
        )


def render_protein_comparison(ref_pdb: str, mut_pdb: str, animate=True):
    """
    Legacy: side-by-side comparison. Kept for compatibility.
    Prefer render_protein_overlay() for the HG-DT wizard UI.
    """
    col1, col2 = st.columns(2)

    with col1:
        st.caption("Reference")
        if ref_pdb:
            show_protein_3d(ref_pdb, width=400, height=380, animate=animate,
                            color=_REF_COLOR)
        else:
            st.info("Reference PDB not provided.")

    with col2:
        st.caption("Mutant")
        if mut_pdb:
            show_protein_3d(mut_pdb, width=400, height=380, animate=animate,
                            color=_MUT_COLOR)
        else:
            st.info("Mutant PDB not provided.")

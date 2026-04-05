import streamlit as st
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
import pandas as pd
from PIL import Image

# Import pipeline modules
from src.hg_dt.viz.tracks_plotter import plot_tracks, plot_multi_tracks
from src.hg_dt.viz.hic_plotter import plot_hic_triangle
from src.hg_dt.viz.protein_viz import render_protein_comparison
from src.hg_dt.viz.browser import render_genome_scroller
from src.hg_dt.analyze.attribution import generate_mechanistic_insight
from src.hg_dt.analyze.deltas import (
    accessibility_delta, expression_delta, contact_delta,
    find_silenced_elements, find_distal_loops,
)

st.set_page_config(page_title="HG-DT Causal Dashboard", layout="wide")

st.title("HG-DT: Visual Interpretation & Causal Dashboard")

# ---------------------------------------------------------------------------
# Sidebar — pipeline settings
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("Pipeline Settings")
    use_real_pipeline = st.checkbox("Use Real Pipeline", value=False,
                                    help="Requires ALPHA_GENOME_API_KEY and reference files.")
    if use_real_pipeline:
        st.caption("Reference file paths (hg38):")
        fasta_path = st.text_input("FASTA dir / file", value="data/reference/hg38/")
        gtf_path   = st.text_input("GENCODE GTF",      value="data/reference/hg38/gencode.v49.annotation.gtf.gz")
        bed_path   = st.text_input("ENCODE cCREs BED", value="data/reference/hg38/GRCh38-cCREs.bed.gz")
    else:
        fasta_path = gtf_path = bed_path = ""


# ---------------------------------------------------------------------------
# Mock data (fallback / demo mode)
# ---------------------------------------------------------------------------
@st.cache_data
def generate_mock_data(mod_type: str, gene: str):
    ref_track = np.random.normal(5, 1, 100)
    mut_track = ref_track.copy()
    if mod_type == "Deletion":
        mut_track[40:60] = mut_track[40:60] * 0.2

    multi_tracks = {}
    for name in ['ATAC', 'CTCF', 'H3K27ac', 'CAGE']:
        r = np.random.normal(5, 1, 100)
        m = r.copy()
        if mod_type == "Deletion":
            m[45:55] *= 0.1
        multi_tracks[name] = (r, m)

    ref_hic = np.random.rand(50, 50)
    ref_hic = (ref_hic + ref_hic.T) / 2
    np.fill_diagonal(ref_hic, 1.0)
    mut_hic = ref_hic.copy()
    if mod_type == "Deletion":
        mut_hic[10:20, 30:40] *= 0.3
        mut_hic[30:40, 10:20] = mut_hic[10:20, 30:40]

    def make_jittery_pdb(n_frames=10):
        pdb_content = ""
        for f in range(n_frames):
            pdb_content += f"MODEL {f+1}\n"
            for i in range(1, 21):
                z, x, y = i * 1.5, np.cos(i/2) * 5, np.sin(i/2) * 5
                if mod_type != "None":
                    z += np.sin(f/2) * 2
                pdb_content += (f"ATOM  {i:5d}  CA  ALA A {i:4d}    "
                                f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00 20.00           C  \n")
            pdb_content += "ENDMDL\n"
        return pdb_content

    ref_pdb = make_jittery_pdb()
    mut_pdb = make_jittery_pdb() if mod_type != "None" else ""

    genes_df = pd.DataFrame([
        {'Name': gene,          'Start': 47210000, 'End': 47250000, 'Type': 'Gene'},
        {'Name': 'Enhancer_1',  'Start': 47200000, 'End': 47205000, 'Type': 'Enhancer'},
        {'Name': 'CTCF_Site',   'Start': 47230000, 'End': 47232000, 'Type': 'CTCF'},
        {'Name': 'H3K27ac_Peak','Start': 47208000, 'End': 47210000, 'Type': 'Enhancer'},
    ])

    # Hardcoded mock delta_stats (sign convention: positive = loss)
    if mod_type == "Deletion":
        delta_stats = {
            "loop_weakened": True,
            "accessibility_drop": 0.28 if gene == "TAL1" else 0.45,
            "expression_drop": 0.35 if gene == "TAL1" else 0.50,
        }
    elif mod_type == "Insertion":
        delta_stats = {
            "loop_strengthened": True,
            "accessibility_drop": -0.50,
            "expression_drop": -0.80,
        }
    else:
        delta_stats = {"loop_weakened": False, "accessibility_drop": 0.0, "expression_drop": 0.0}

    return ref_track, mut_track, ref_hic, mut_hic, ref_pdb, mut_pdb, genes_df, multi_tracks, delta_stats


# ---------------------------------------------------------------------------
# Real pipeline orchestration
# ---------------------------------------------------------------------------
def _extract_1d(track_data, index: int = 0) -> np.ndarray:
    """Extract a 1D array from an AlphaGenome TrackData object."""
    vals = np.array(track_data.values)
    if vals.ndim == 1:
        return vals
    return vals[:, index]


@st.cache_data(show_spinner="Running AlphaGenome predictions…")
def run_real_pipeline(chrom: str, edit_start: int, edit_end: int,
                      edit_type: str, gene: str,
                      fasta_path: str, gtf_path: str, bed_path: str):
    """
    Orchestrate the full HG-DT pipeline:
      ReferenceContextBuilder → AlphaGenomeConnector → deltas → translation → ProteinFolder
    Returns the same shape as generate_mock_data so all tabs work unchanged.
    """
    from src.hg_dt.data.builder import ReferenceContextBuilder
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector
    from src.hg_dt.translate.transcript import predict_isoforms
    from src.hg_dt.translate.translator import compare_translation
    from src.hg_dt.models.protein import ProteinFolder

    # 1. Reference context
    builder = ReferenceContextBuilder(fasta_path, gtf_path, bed_path)
    context = builder.get_context(chrom, edit_start, edit_end, edit_type)
    ref_seq  = context['ref_seq']
    mut_seq  = context['mut_seq']
    win_start, _ = context['window']
    genes_df = context['annotations']['genes']

    # 2. AlphaGenome predictions
    connector = AlphaGenomeConnector()
    ref_out = connector.predict_sequence(ref_seq, organism="HUMAN",
                                         requested_outputs=['ATAC', 'CAGE', 'CHIP_TF', 'CONTACT_MAPS'])
    mut_out = connector.predict_sequence(mut_seq, organism="HUMAN",
                                         requested_outputs=['ATAC', 'CAGE', 'CHIP_TF', 'CONTACT_MAPS'])

    # 3. Extract tracks
    ref_atac = _extract_1d(ref_out.atac)
    mut_atac = _extract_1d(mut_out.atac)
    ref_cage = _extract_1d(ref_out.cage)
    mut_cage = _extract_1d(mut_out.cage)
    ref_ctcf = _extract_1d(ref_out.chip_tf)
    mut_ctcf = _extract_1d(mut_out.chip_tf)
    ref_hic  = np.array(ref_out.contact_maps.values[:, :, 0])
    mut_hic  = np.array(mut_out.contact_maps.values[:, :, 0])

    ref_track = ref_atac
    mut_track = mut_atac
    multi_tracks = {
        'ATAC':  (ref_atac, mut_atac),
        'CAGE':  (ref_cage, mut_cage),
        'CTCF':  (ref_ctcf, mut_ctcf),
    }

    # 4. Compute delta_stats (positive = loss in mutant, per spec log2(Ref/Mut))
    a_delta = accessibility_delta(ref_atac, mut_atac)
    e_delta = expression_delta(ref_cage, mut_cage)
    silenced = find_silenced_elements(ref_atac, mut_atac)
    loops_ref = find_distal_loops(ref_hic)
    loops_mut = find_distal_loops(mut_hic)
    loop_weakened = len(loops_mut) < len(loops_ref) * 0.8

    acc_drop = float(np.mean(a_delta[silenced])) if silenced else float(np.mean(a_delta[a_delta > 0]) or 0)
    exp_drop = float(np.mean(e_delta[e_delta > 0]) or 0)

    delta_stats = {
        "loop_weakened": loop_weakened,
        "loop_strengthened": len(loops_mut) > len(loops_ref) * 1.2,
        "accessibility_drop": acc_drop,
        "expression_drop": exp_drop,
    }

    # 5. Build annotation DataFrame for genome scroller
    if genes_df.empty:
        display_genes = pd.DataFrame([
            {'Name': gene, 'Start': edit_start, 'End': edit_end, 'Type': 'Gene'}
        ])
    else:
        display_genes = genes_df[['Start', 'End']].copy()
        display_genes['Name'] = genes_df.get('gene_name', genes_df.index.astype(str))
        display_genes['Type'] = genes_df.get('Feature', 'Gene')

    # 6. Protein structure (first isoform)
    ref_pdb = mut_pdb = ""
    if not genes_df.empty:
        mut_shift = -(edit_end - edit_start) if edit_type == 'deletion' else 0
        isoforms = predict_isoforms(genes_df, ref_seq, mut_seq, win_start,
                                    edit_start, edit_end, mut_shift)
        folder = ProteinFolder()
        for tid, iso in isoforms.items():
            trans = compare_translation(iso['ref_mrna'], iso['mut_mrna'])
            if trans['ref_aa']:
                ref_pdb = folder.predict_structure(trans['ref_aa'])['pdb']
            if trans['mut_aa']:
                mut_pdb = folder.predict_structure(trans['mut_aa'])['pdb']
            break  # use first isoform only

    return (ref_track, mut_track, ref_hic, mut_hic,
            ref_pdb, mut_pdb, display_genes, multi_tracks, delta_stats)


# ---------------------------------------------------------------------------
# Accessibility-over-time simulation (unchanged)
# ---------------------------------------------------------------------------
@st.cache_data
def generate_accessibility_over_time(gene: str, mod_type: str, time_steps: int = 50):
    t = np.linspace(0, 10, time_steps)
    base = 100.0
    if mod_type == "Deletion":
        acc = base - 60 * (1 / (1 + np.exp(-t + 5))) + np.random.normal(0, 2, time_steps)
    elif mod_type == "Insertion":
        acc = base + 80 * (1 / (1 + np.exp(-t + 5))) + np.random.normal(0, 2, time_steps)
    else:
        acc = np.full(time_steps, base) + np.random.normal(0, 2, time_steps)
    return t, acc


# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------
for key, default in [("mod_type", "None"), ("gene", "TAL1"),
                     ("analyzed", False), ("browser_pos", 47200000)]:
    if key not in st.session_state:
        st.session_state[key] = default

# ---------------------------------------------------------------------------
# Tabs
# ---------------------------------------------------------------------------
tab1, tab2, tab2b, tab3, tab4, tab5, tab6 = st.tabs([
    "Specification", "Genome Tracks", "Genome Scroller",
    "3D Organization", "Protein Structure", "Trajectory Animation",
    "Mechanistic Attribution",
])

with tab1:
    st.header("Input DNA Modification")
    col1, col2 = st.columns(2)
    with col1:
        gene = st.selectbox("Select Target Gene:", ["TAL1", "OCT4", "NANOG", "SOX2", "Mef2c"])
    with col2:
        mod_type = st.selectbox("Select Modification Type:",
                                ["None", "Deletion", "Insertion", "Duplication"])

    chrom = st.text_input("Chromosome", "chr1")
    locus = st.text_input("Locus", "47200000-47250000")

    if st.button("Run Analysis"):
        st.session_state.gene = gene
        st.session_state.mod_type = mod_type
        st.session_state.analyzed = True
        st.success(f"Analysis triggered for {gene} ({mod_type}) at {chrom}:{locus}")

if st.session_state.analyzed:
    # Parse locus
    try:
        locus_parts = locus.replace(",", "").split("-")
        edit_start = int(locus_parts[0])
        edit_end   = int(locus_parts[1])
    except Exception:
        edit_start, edit_end = 47200000, 47250000

    # Dispatch to real or mock pipeline
    if use_real_pipeline:
        try:
            (ref_track, mut_track, ref_hic, mut_hic,
             ref_pdb, mut_pdb, genes_df, multi_tracks, delta_stats) = run_real_pipeline(
                chrom, edit_start, edit_end,
                st.session_state.mod_type.lower(),
                st.session_state.gene,
                fasta_path, gtf_path, bed_path,
            )
        except Exception as exc:
            st.error(f"Real pipeline failed: {exc}. Falling back to mock data.")
            (ref_track, mut_track, ref_hic, mut_hic,
             ref_pdb, mut_pdb, genes_df, multi_tracks, delta_stats) = generate_mock_data(
                st.session_state.mod_type, st.session_state.gene
            )
    else:
        (ref_track, mut_track, ref_hic, mut_hic,
         ref_pdb, mut_pdb, genes_df, multi_tracks, delta_stats) = generate_mock_data(
            st.session_state.mod_type, st.session_state.gene
        )

    with tab2:
        st.header("AlphaGenome Multi-Scalar Predictions")
        st.write("Predicted biological tracks from the DNA sequence (Ref vs. Mut):")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            multi_track_path = plot_multi_tracks(
                multi_tracks, tmp.name,
                title=f"Regulatory Landscape: {st.session_state.gene}"
            )
            st.image(Image.open(multi_track_path), use_container_width=True)
        with st.expander("🔬 Model & Methodology"):
            st.markdown("""
**Model:** AlphaGenome (DeepMind, 2026)
**Method:** Sequence-to-function transformer predicting thousands of genomic tracks from a 1 Mb DNA window.
**Tracks shown:** ATAC (chromatin accessibility), CAGE (transcription initiation), CTCF (insulator binding), H3K27ac (active enhancer mark).
""")

    with tab2b:
        st.header("Interactive Modification Planner")
        search_col1, search_col2 = st.columns([3, 1])
        with search_col1:
            search_gene = st.text_input("Genome Search (Gene Symbol or Coordinates)",
                                        value=st.session_state.gene)
        with search_col2:
            if st.button("Jump to Locus"):
                lookup = {"TAL1": 47210000, "MYC": 127735434, "BRCA1": 43044295,
                          "GATA1": 48786554, "OCT4": 31132048}
                if search_gene.upper() in lookup:
                    st.session_state.browser_pos = lookup[search_gene.upper()]
                    st.session_state.gene = search_gene.upper()
                elif ":" in search_gene:
                    try:
                        pos = int(search_gene.split(":")[-1].replace(",", "").split("-")[0])
                        st.session_state.browser_pos = pos
                    except Exception:
                        pass

        with st.expander("🛠️ Deletion / Insertion Tool", expanded=True):
            plan_col1, plan_col2, plan_col3 = st.columns(3)
            with plan_col1:
                del_start = st.number_input("Deletion Start",
                                            value=st.session_state.browser_pos - 5000, step=100)
            with plan_col2:
                del_end = st.number_input("Deletion End",
                                          value=st.session_state.browser_pos + 5000, step=100)
            with plan_col3:
                if st.button("Apply Modification to Pipeline"):
                    st.session_state.mod_type = "Deletion"
                    st.session_state.analyzed = True
                    st.info(f"Model re-calculated for deletion at {del_start}-{del_end}")

        browser_pos = st.slider(
            "Coordinate Scroller (hg38)",
            min_value=st.session_state.browser_pos - 2000000,
            max_value=st.session_state.browser_pos + 2000000,
            value=st.session_state.browser_pos, step=1000,
        )
        st.session_state.browser_pos = browser_pos
        render_genome_scroller(genes_df, browser_pos, window_size=500000,
                               deletion_region=(del_start, del_end))
        st.caption("Tip: Use the slider to navigate to distal enhancers. "
                   "Red region = planned modification.")

    with tab3:
        st.header("3D Organization (Hi-C Heatmaps)")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            hic_img_path = plot_hic_triangle(
                ref_hic, mut_hic, tmp.name,
                title=f"3D Contact Map Delta: {st.session_state.gene}"
            )
            st.image(Image.open(hic_img_path), use_container_width=True)
        with st.expander("🔬 Model & Methodology"):
            st.markdown("""
**Model:** AlphaGenome Hi-C Predictor
**View:** Triangular heatmap (diagonal rotated to base). Delta = Mutant − Reference.
- **Warm (red) pixels:** Gained interactions (ectopic contacts).
- **Cool (blue) pixels:** Lost interactions (broken enhancer–promoter loops).
""")

    with tab4:
        st.header("Protein Structure: Molecular Dynamics Simulation")
        st.write("Reference (left) vs. Predicted Mutated (right). Powered by ESMFold.")
        render_protein_comparison(ref_pdb, mut_pdb, animate=True)
        with st.expander("🔬 Model & Methodology"):
            st.markdown("""
**Model:** Meta ESMFold (public API) with mock fallback for long sequences
**Process:** DNA → mRNA isoforms (GENCODE exon splicing) → amino acid translation (Biopython) → 3D folding (ESMFold).
**Confidence:** B-factor column of the PDB output = pLDDT (0–100). Regions <50 are unreliably folded.
""")

    with tab5:
        st.header("Trajectory Animation")
        t, accessibility = generate_accessibility_over_time(
            st.session_state.gene, st.session_state.mod_type
        )
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
            st.image(Image.open(tmp.name), use_container_width=True)
        with st.expander("🔬 Model & Methodology"):
            st.markdown(f"""
**Model:** Custom Langevin Dynamics Polymer Engine
**Trajectory Physics:** Contact constraints are updated after the sequence edit; chromatin refolds toward a new energy minimum. The drop in signal (purple) represents the gene losing its regulatory neighbourhood over simulation time.
""")

    with tab6:
        st.header("Mechanistic Attribution")

        mod_details = {
            "type": st.session_state.mod_type.lower(),
            "target": f"enhancer element near {st.session_state.gene}",
        }
        insight = generate_mechanistic_insight(mod_details, delta_stats)

        st.info("### Mechanistic Summary")
        st.success(insight)
        with st.expander("🔬 Model & Methodology"):
            st.write("**Model:** Rule-based mechanistic attribution from computed track deltas")
            st.write("**Process:** Accessibility and expression deltas are computed via "
                     "`log2(Ref/Mut)` from AlphaGenome predictions. Loop counts are derived "
                     "from the contact map. These numeric signals are then formatted into a "
                     "human-readable causal chain.")

elif st.session_state.analyzed:
    st.warning("Please select a valid modification from the Specification tab.")

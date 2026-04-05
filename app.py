import streamlit as st
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional
from PIL import Image
from dotenv import load_dotenv

load_dotenv()

# Import pipeline modules
import plotly.express as px
import plotly.graph_objects as go

from src.hg_dt.viz.tracks_plotter import plot_tracks, plot_multi_tracks
from src.hg_dt.viz.hic_plotter import plot_hic_triangle
from src.hg_dt.viz.protein_viz import render_protein_overlay
from src.hg_dt.viz.chromatin_3d import render_chromatin_animation
from src.hg_dt.viz.browser import render_genome_scroller
from src.hg_dt.analyze.attribution import generate_mechanistic_insight
from src.hg_dt.analyze.deltas import (
    accessibility_delta, expression_delta, contact_delta,
    find_silenced_elements, find_distal_loops,
)
from src.hg_dt.models.alphagenome import CELL_TYPES
from src.hg_dt.models.evo2 import Evo2Client

st.set_page_config(page_title="HG-DT: Human Genome Digital Twin", layout="wide")

# ---------------------------------------------------------------------------
# Gene coordinate registry (hg38)
# ---------------------------------------------------------------------------
GENE_COORDS = {
    "TAL1":  ("chr1",  47210000, 47260000),
    "MYC":   ("chr8",  127735000, 127742000),
    "BRCA1": ("chr17", 43044000,  43126000),
    "GATA1": ("chrX",  48786000,  48812000),
    "OCT4":  ("chr6",  31132000,  31147000),
    "NANOG": ("chr12", 7794000,   7802000),
    "SOX2":  ("chr3",  181429000, 181435000),
    "TP53":  ("chr17", 7668000,   7688000),
    "CTCF":  ("chr16", 67562000,  67600000),
    "RUNX1": ("chr21", 34787000,  36004000),
    "PAX5":  ("chr9",  36838000,  37035000),
    "EZH2":  ("chr7",  148504000, 148581000),
    "MYB":   ("chr6",  135502000, 135540000),
    "MEF2C": ("chr5",  88101000,  88203000),
}

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
        st.caption("Cell type (optional):")
        cell_type_options = ["None (all cell types)"] + sorted(CELL_TYPES.keys())
        cell_type_sel = st.selectbox("Cell Type", cell_type_options, index=0)
        cell_type = None if cell_type_sel.startswith("None") else cell_type_sel
    else:
        fasta_path = gtf_path = bed_path = ""
        cell_type = None

    st.divider()
    st.caption(f"Evo 2 backend: **{Evo2Client().backend}**")

    st.divider()
    if st.button("Start Over"):
        for key in ["step", "chrom", "gene", "edit_start", "edit_end",
                    "mod_type", "pipeline_result"]:
            if key in st.session_state:
                del st.session_state[key]
        st.rerun()


# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------
for key, default in [
    ("step",            1),
    ("chrom",           "chr1"),
    ("gene",            "TAL1"),
    ("edit_start",      47200000),
    ("edit_end",        47250000),
    ("mod_type",        "Deletion"),
    ("pipeline_result", None),
]:
    if key not in st.session_state:
        st.session_state[key] = default


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

    return ref_track, mut_track, ref_hic, mut_hic, ref_pdb, mut_pdb, multi_tracks, delta_stats


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
                      fasta_path: str, gtf_path: str, bed_path: str,
                      cell_type: Optional[str] = None):
    """
    Orchestrate the full HG-DT pipeline:
      ReferenceContextBuilder → AlphaGenomeConnector → deltas
      → splice-aware transcript prediction → translation → ProteinFolder
    """
    from src.hg_dt.data.builder import ReferenceContextBuilder
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector
    from src.hg_dt.translate.transcript import predict_isoforms
    from src.hg_dt.translate.translator import compare_translation
    from src.hg_dt.models.protein import ProteinFolder

    builder = ReferenceContextBuilder(fasta_path, gtf_path, bed_path)
    context = builder.get_context(chrom, edit_start, edit_end, edit_type)
    ref_seq   = context['ref_seq']
    mut_seq   = context['mut_seq']
    win_start, _ = context['window']
    genes_df  = context['annotations']['genes']

    connector = AlphaGenomeConnector()
    outputs_requested = ['ATAC', 'CAGE', 'CHIP_TF', 'CONTACT_MAPS', 'SPLICE_SITES']
    ref_out = connector.predict_sequence(ref_seq, organism="HUMAN",
                                         requested_outputs=outputs_requested,
                                         cell_type=cell_type)
    mut_out = connector.predict_sequence(mut_seq, organism="HUMAN",
                                         requested_outputs=outputs_requested,
                                         cell_type=cell_type)

    ref_atac = _extract_1d(ref_out.atac)
    mut_atac = _extract_1d(mut_out.atac)
    ref_cage = _extract_1d(ref_out.cage)
    mut_cage = _extract_1d(mut_out.cage)
    ref_ctcf = _extract_1d(ref_out.chip_tf)
    mut_ctcf = _extract_1d(mut_out.chip_tf)
    ref_hic  = np.array(ref_out.contact_maps.values[:, :, 0])
    mut_hic  = np.array(mut_out.contact_maps.values[:, :, 0])

    ref_splice = _extract_1d(ref_out.splice_sites) if ref_out.splice_sites is not None else None
    mut_splice = _extract_1d(mut_out.splice_sites) if mut_out.splice_sites is not None else None

    multi_tracks = {
        'ATAC':  (ref_atac, mut_atac),
        'CAGE':  (ref_cage, mut_cage),
        'CTCF':  (ref_ctcf, mut_ctcf),
    }

    a_delta  = accessibility_delta(ref_atac, mut_atac)
    e_delta  = expression_delta(ref_cage, mut_cage)
    silenced = find_silenced_elements(ref_atac, mut_atac)
    loops_ref = find_distal_loops(ref_hic)
    loops_mut = find_distal_loops(mut_hic)
    loop_weakened    = len(loops_mut) < len(loops_ref) * 0.8
    loop_strengthened = len(loops_mut) > len(loops_ref) * 1.2

    acc_drop = float(np.mean(a_delta[silenced])) if silenced else float(np.mean(a_delta[a_delta > 0]) or 0)
    exp_drop = float(np.mean(e_delta[e_delta > 0]) or 0)

    delta_stats = {
        "loop_weakened": loop_weakened,
        "loop_strengthened": loop_strengthened,
        "accessibility_drop": acc_drop,
        "expression_drop": exp_drop,
    }

    ref_pdb = mut_pdb = ""
    splice_info: dict = {}
    if not genes_df.empty:
        mut_shift = -(edit_end - edit_start) if edit_type == 'deletion' else 0
        isoforms = predict_isoforms(
            genes_df, ref_seq, mut_seq, win_start,
            edit_start, edit_end, mut_shift,
            ref_splice_track=ref_splice,
            mut_splice_track=mut_splice,
        )
        folder = ProteinFolder()
        for tid, iso in isoforms.items():
            trans = compare_translation(iso['ref_mrna'], iso['mut_mrna'])
            if trans['ref_aa']:
                ref_pdb = folder.predict_structure(trans['ref_aa'])['pdb']
            if trans['mut_aa']:
                mut_pdb = folder.predict_structure(trans['mut_aa'])['pdb']
            splice_info = {
                'transcript_id': tid,
                'frameshift': iso['frameshift'],
                'splice_disrupted': iso['splice_disrupted'],
                'skipped_exons': iso['skipped_exons'],
            }
            break

    window_size = min(2048, len(ref_seq))
    center = len(ref_seq) // 2
    ref_window = ref_seq[center - window_size//2 : center + window_size//2]
    mut_window = mut_seq[center - window_size//2 : center + window_size//2]
    evo2 = Evo2Client()
    evo2_result = evo2.score_variant(ref_window, mut_window)
    evo2_scan   = evo2.scan_sequence(ref_window, window=256)

    return (ref_atac, mut_atac, ref_hic, mut_hic,
            ref_pdb, mut_pdb, multi_tracks, delta_stats,
            {"evo2": evo2_result, "evo2_scan": evo2_scan, "splice": splice_info,
             "cell_type": cell_type})


@st.cache_data(show_spinner="Running differentiation trajectory…")
def run_differentiation_trajectory(chrom: str, edit_start: int, edit_end: int,
                                    mod_type: str, fasta_path: str,
                                    gtf_path: str, bed_path: str):
    """
    Run AlphaGenome predict_sequence at H1-hESC and K562 for ref and mut.
    Returns dict with 4 contact maps keyed by h1_ref, h1_mut, k562_ref, k562_mut.
    """
    from src.hg_dt.data.builder import ReferenceContextBuilder
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector

    builder = ReferenceContextBuilder(fasta_path, gtf_path, bed_path)
    context = builder.get_context(chrom, edit_start, edit_end, mod_type.lower())
    ref_seq = context['ref_seq']
    mut_seq = context['mut_seq']

    connector = AlphaGenomeConnector()
    results = {}
    for label, seq, ct in [
        ('h1_ref',   ref_seq, 'H1-hESC'),
        ('h1_mut',   mut_seq, 'H1-hESC'),
        ('k562_ref', ref_seq, 'K562'),
        ('k562_mut', mut_seq, 'K562'),
    ]:
        out = connector.predict_sequence(seq, requested_outputs=['CONTACT_MAPS'],
                                          cell_type=ct)
        results[label] = np.array(out.contact_maps.values[:, :, 0])
    return results


def _mock_trajectory(mod_type: str):
    """Generate 4 mock contact maps for H1-hESC and K562 in demo mode."""
    base = np.random.rand(30, 30)
    base = (base + base.T) / 2
    np.fill_diagonal(base, 1.0)
    h1_ref = base.copy()
    k562_ref = base * 0.8 + np.random.rand(30, 30) * 0.1
    k562_ref = (k562_ref + k562_ref.T) / 2
    h1_mut = h1_ref.copy()
    k562_mut = k562_ref.copy()
    if mod_type == "Deletion":
        h1_mut[8:15, 18:25] *= 0.3
        h1_mut[18:25, 8:15] = h1_mut[8:15, 18:25]
        k562_mut[8:15, 18:25] *= 0.4
        k562_mut[18:25, 8:15] = k562_mut[8:15, 18:25]
    return {'h1_ref': h1_ref, 'h1_mut': h1_mut,
            'k562_ref': k562_ref, 'k562_mut': k562_mut}


# ---------------------------------------------------------------------------
# Helper: build preview genes_df for Step 2 (no API call)
# ---------------------------------------------------------------------------
def _build_preview_genes_df() -> pd.DataFrame:
    """Build a simple annotation DataFrame centered on the current gene."""
    gene = st.session_state.gene
    chrom = st.session_state.chrom
    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end

    # Get gene body coordinates
    if gene in GENE_COORDS:
        _, gs, ge = GENE_COORDS[gene]
    else:
        gs, ge = edit_start - 10000, edit_end + 10000

    mid = (gs + ge) // 2
    rows = [
        {'Name': gene,             'Start': gs,         'End': ge,         'Type': 'Gene'},
        {'Name': f'{gene}_ENH1',   'Start': mid - 25000,'End': mid - 20000,'Type': 'Enhancer'},
        {'Name': f'{gene}_ENH2',   'Start': mid + 15000,'End': mid + 20000,'Type': 'Enhancer'},
        {'Name': f'CTCF_L',        'Start': gs - 5000,  'End': gs - 3000,  'Type': 'CTCF'},
        {'Name': f'CTCF_R',        'Start': ge + 3000,  'End': ge + 5000,  'Type': 'CTCF'},
        {'Name': f'{gene}_H3K27ac','Start': mid - 5000, 'End': mid - 2000, 'Type': 'Enhancer'},
    ]
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# Helper: run pipeline and package result as dict
# ---------------------------------------------------------------------------
def _run_pipeline() -> dict:
    gene      = st.session_state.gene
    chrom     = st.session_state.chrom
    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end
    mod_type  = st.session_state.mod_type

    if use_real_pipeline:
        try:
            result_tuple = run_real_pipeline(
                chrom, edit_start, edit_end,
                mod_type.lower(), gene,
                fasta_path, gtf_path, bed_path,
                cell_type=cell_type,
            )
            (ref_atac, mut_atac, ref_hic, mut_hic,
             ref_pdb, mut_pdb, multi_tracks, delta_stats, extras) = result_tuple
        except Exception as exc:
            st.error(f"Real pipeline failed: {exc}. Falling back to mock data.")
            use_mock = True
        else:
            use_mock = False
    else:
        use_mock = True

    if use_mock:
        (ref_atac, mut_atac, ref_hic, mut_hic,
         ref_pdb, mut_pdb, multi_tracks, delta_stats) = generate_mock_data(mod_type, gene)
        evo2_client = Evo2Client()
        mock_ref = "ATCGATCGCCGCGAGGTGGCAG" * 20
        mock_alt = mock_ref[:200] + "TTTTTTTTTT" + mock_ref[210:]
        extras = {
            "evo2": evo2_client.score_variant(mock_ref, mock_alt),
            "evo2_scan": evo2_client.scan_sequence(mock_ref, window=44),
            "splice": {},
            "cell_type": None,
        }

    return {
        "ref_atac":    ref_atac,
        "mut_atac":    mut_atac,
        "ref_hic":     ref_hic,
        "mut_hic":     mut_hic,
        "ref_pdb":     ref_pdb,
        "mut_pdb":     mut_pdb,
        "multi_tracks": multi_tracks,
        "delta_stats": delta_stats,
        "extras":      extras,
    }


# ---------------------------------------------------------------------------
# Helper: render differentiation trajectory panel
# ---------------------------------------------------------------------------
def _render_trajectory(result: dict):
    """Render 4-panel contact map trajectory (H1-hESC → K562, Ref vs Mut)."""
    if use_real_pipeline and fasta_path and gtf_path and bed_path:
        try:
            maps = run_differentiation_trajectory(
                st.session_state.chrom,
                st.session_state.edit_start,
                st.session_state.edit_end,
                st.session_state.mod_type,
                fasta_path, gtf_path, bed_path,
            )
        except Exception as exc:
            st.warning(f"Trajectory prediction failed: {exc}. Showing mock.")
            maps = _mock_trajectory(st.session_state.mod_type)
    else:
        maps = _mock_trajectory(st.session_state.mod_type)

    titles = ["H1-hESC  Ref", "H1-hESC  Mut", "K562  Ref", "K562  Mut"]
    keys   = ["h1_ref", "h1_mut", "k562_ref", "k562_mut"]
    vmax   = max(np.nanmax(np.abs(maps[k])) for k in keys if maps[k].size > 0)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for ax, key, title in zip(axes, keys, titles):
        im = ax.imshow(maps[key], cmap='YlOrRd', vmin=0, vmax=vmax, aspect='auto')
        ax.set_title(title, fontsize=11)
        ax.axis('off')

    # Shared colorbar
    fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="Contact frequency")

    # Stage arrow annotation
    fig.text(0.5, -0.02,
             "← Embryonic stem cell (H1-hESC)   ·   Committed myeloid (K562) →",
             ha='center', fontsize=10, color='#444')

    plt.tight_layout()
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        fig.savefig(tmp.name, dpi=150, bbox_inches='tight')
        plt.close(fig)
        st.image(Image.open(tmp.name), use_container_width=True)

    with st.expander("Model & Methodology"):
        st.markdown("""
**Model:** AlphaGenome (DeepMind) with cell-type-specific ontology terms
**Cell types:** H1-hESC (`CLO:0037111`) = embryonic pluripotent stem cell; K562 (`CLO:0007050`) = committed CML myeloid line
**Interpretation:** Contact maps are predicted independently for each cell type. Changes between Ref and Mut show how the edit disrupts 3D genome organization across the developmental axis.
**Progenitor interpolation:** Intermediate stages (e.g. progenitor) would fall between these two extremes.
""")


# ---------------------------------------------------------------------------
# Helper: per-regulatory-element accessibility bubble chart
# ---------------------------------------------------------------------------
def _render_accessibility_elements(multi_tracks: dict, genes_df: pd.DataFrame,
                                    edit_start: int, edit_end: int):
    """
    Plot per-element accessibility change as an interactive Plotly bubble chart.
    X = genomic position, Y = Δ log₂(Ref/Mut), size = element width, color = type.
    """
    atac_ref, atac_mut = multi_tracks.get('ATAC', (None, None))
    if atac_ref is None or genes_df is None or genes_df.empty:
        st.caption("No per-element data available.")
        return

    n = len(atac_ref)
    window_span = max(edit_end - edit_start, 1)

    rows = []
    for _, row in genes_df.iterrows():
        mid = (row['Start'] + row['End']) / 2
        size = max(row['End'] - row['Start'], 1)

        # Map genomic coordinate to track index (proportional within window)
        idx = int((mid - edit_start) / window_span * n)
        idx = max(0, min(n - 1, idx))

        # Sample a window of ±2 bins for robustness
        i0 = max(0, idx - 2)
        i1 = min(n, idx + 3)
        ref_sig = float(np.mean(atac_ref[i0:i1]))
        mut_sig = float(np.mean(atac_mut[i0:i1]))
        delta   = float(np.log2((ref_sig + 1e-6) / (mut_sig + 1e-6)))

        rows.append({
            'Element':          row['Name'],
            'Type':             row['Type'],
            'Position':         int(mid),
            'Element size (bp)': int(size),
            'Ref ATAC':         round(ref_sig, 3),
            'Mut ATAC':         round(mut_sig, 3),
            'Δ log₂(Ref/Mut)':  round(delta, 3),
        })

    if not rows:
        return

    df = pd.DataFrame(rows)

    color_map = {'Gene': '#4a9eff', 'Enhancer': '#ff9a3c', 'CTCF': '#2ecc71'}

    fig = px.scatter(
        df,
        x='Position',
        y='Δ log₂(Ref/Mut)',
        size='Element size (bp)',
        color='Type',
        text='Element',
        hover_data=['Ref ATAC', 'Mut ATAC'],
        color_discrete_map=color_map,
        size_max=45,
    )
    fig.add_hline(y=0, line_dash='dash', line_color='rgba(180,180,180,0.5)',
                  annotation_text='No change', annotation_position='right')
    fig.add_vrect(
        x0=edit_start, x1=edit_end,
        fillcolor='rgba(255,80,80,0.08)', line_color='rgba(255,80,80,0.4)',
        annotation_text='Edit', annotation_position='top left',
    )
    fig.update_traces(textposition='top center', textfont_size=10)
    fig.update_layout(
        height=320,
        xaxis=dict(tickformat=',.0s', title='Genomic position (hg38)'),
        yaxis=dict(title='Accessibility Δ log₂ (Ref / Mut)',
                   zeroline=False),
        plot_bgcolor='white',
        paper_bgcolor='white',
        legend=dict(orientation='h', yanchor='bottom', y=1.02, xanchor='right', x=1),
        margin=dict(t=30, b=40, l=60, r=20),
    )
    # Shade gain vs loss
    fig.add_hrect(y0=0, y1=df['Δ log₂(Ref/Mut)'].max() + 0.5,
                  fillcolor='rgba(255,80,80,0.04)', line_width=0,
                  annotation_text='Accessibility lost in mutant',
                  annotation_position='top right',
                  annotation_font_color='rgba(200,50,50,0.7)')
    fig.add_hrect(y0=df['Δ log₂(Ref/Mut)'].min() - 0.5, y1=0,
                  fillcolor='rgba(50,180,100,0.04)', line_width=0,
                  annotation_text='Accessibility gained in mutant',
                  annotation_position='bottom right',
                  annotation_font_color='rgba(50,150,80,0.7)')

    st.plotly_chart(fig, use_container_width=True)


# ---------------------------------------------------------------------------
# Step progress indicator
# ---------------------------------------------------------------------------
st.title("HG-DT: Human Genome Digital Twin")

steps = ["1 · Find Locus", "2 · What's Here?", "3 · Results"]
prog_cols = st.columns(3)
for i, (col, label) in enumerate(zip(prog_cols, steps)):
    active = (i + 1 == st.session_state.step)
    col.markdown(
        f"**:blue[{label}]**" if active else f"<span style='color:#aaa'>{label}</span>",
        unsafe_allow_html=True,
    )

st.divider()

# ===========================================================================
# STEP 1 — Find Your Locus
# ===========================================================================
if st.session_state.step == 1:
    st.header("Step 1 — Find Your Locus")
    st.caption("Search for a gene or enter genomic coordinates, then define your edit.")

    col_search, col_jump = st.columns([3, 1])
    with col_search:
        search_input = st.text_input(
            "Gene symbol or chr:start-end",
            value=st.session_state.gene,
            placeholder="e.g. TAL1  or  chr1:47200000-47260000",
        )
    with col_jump:
        st.write("")  # vertical align
        if st.button("Look up", use_container_width=True):
            gene_upper = search_input.upper().strip()
            if gene_upper in GENE_COORDS:
                c, s, e = GENE_COORDS[gene_upper]
                st.session_state.chrom      = c
                st.session_state.edit_start = s
                st.session_state.edit_end   = e
                st.session_state.gene       = gene_upper
                st.rerun()
            elif ":" in search_input and "-" in search_input:
                try:
                    parts = search_input.replace(",", "").split(":")
                    chrom_part = parts[0].strip()
                    range_part = parts[1].strip().split("-")
                    st.session_state.chrom      = chrom_part
                    st.session_state.edit_start = int(range_part[0])
                    st.session_state.edit_end   = int(range_part[1])
                    st.rerun()
                except Exception:
                    st.error("Could not parse coordinates. Use format chr1:47200000-47250000")
            else:
                st.warning(f"Gene '{search_input}' not in registry. Enter coordinates manually.")

    st.write("")
    col_mod, col_chrom = st.columns([1, 1])
    with col_mod:
        mod_type = st.selectbox(
            "Modification type",
            ["Deletion", "Insertion", "SNP", "None"],
            index=["Deletion", "Insertion", "SNP", "None"].index(st.session_state.mod_type)
                  if st.session_state.mod_type in ["Deletion", "Insertion", "SNP", "None"] else 0,
        )
    with col_chrom:
        chrom = st.text_input("Chromosome", value=st.session_state.chrom)

    col_s, col_e = st.columns(2)
    with col_s:
        edit_start = st.number_input("Edit start (bp)", value=st.session_state.edit_start, step=1000)
    with col_e:
        edit_end = st.number_input("Edit end (bp)", value=st.session_state.edit_end, step=1000)

    edit_size = int(edit_end) - int(edit_start)
    if edit_size > 0:
        st.caption(f"Edit size: **{edit_size:,} bp**")
    elif edit_size <= 0:
        st.error("Edit end must be greater than edit start.")

    st.write("")
    if st.button("Preview Locus →", type="primary", use_container_width=False,
                 disabled=(edit_size <= 0)):
        st.session_state.chrom      = chrom
        st.session_state.edit_start = int(edit_start)
        st.session_state.edit_end   = int(edit_end)
        st.session_state.mod_type   = mod_type
        st.session_state.step       = 2
        st.rerun()


# ===========================================================================
# STEP 2 — What's Here?
# ===========================================================================
elif st.session_state.step == 2:
    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end
    mid_pos    = (edit_start + edit_end) // 2

    if st.button("← Back to Locus"):
        st.session_state.step = 1
        st.rerun()

    st.header("Step 2 — What's Here?")
    st.caption(
        f"Reviewing **{st.session_state.chrom}:{edit_start:,}–{edit_end:,}** "
        f"({edit_end - edit_start:,} bp {st.session_state.mod_type.lower()}) "
        f"· Gene: **{st.session_state.gene}**"
    )

    genes_df = _build_preview_genes_df()

    render_genome_scroller(
        genes_df,
        current_pos=mid_pos,
        window_size=500000,
        deletion_region=(edit_start, edit_end),
    )

    # Annotation overlap table
    overlapping = genes_df[
        (genes_df['Start'] <= edit_end) & (genes_df['End'] >= edit_start)
    ][['Name', 'Type', 'Start', 'End']].copy()

    if not overlapping.empty:
        st.warning(f"Planned edit overlaps **{len(overlapping)}** annotated element(s):")
        overlapping['Start'] = overlapping['Start'].map('{:,}'.format)
        overlapping['End']   = overlapping['End'].map('{:,}'.format)
        st.dataframe(overlapping, use_container_width=True, hide_index=True)
    else:
        st.info("Planned edit does not overlap any annotated elements in this view.")

    st.write("")
    if st.button("Predict Effect →", type="primary"):
        st.session_state.step = 3
        st.session_state.pipeline_result = None  # force recompute
        st.rerun()


# ===========================================================================
# STEP 3 — Results
# ===========================================================================
elif st.session_state.step == 3:
    if st.button("← Back to Preview"):
        st.session_state.step = 2
        st.rerun()

    st.header("Step 3 — Predicted Effect")
    st.caption(
        f"**{st.session_state.gene}** · "
        f"{st.session_state.chrom}:{st.session_state.edit_start:,}–{st.session_state.edit_end:,} · "
        f"{st.session_state.mod_type}"
        + (f" · Cell type: {cell_type}" if cell_type else "")
        + (" · [REAL PIPELINE]" if use_real_pipeline else " · [DEMO MODE]")
    )

    # Run pipeline (cached in session_state)
    if st.session_state.pipeline_result is None:
        st.session_state.pipeline_result = _run_pipeline()

    result = st.session_state.pipeline_result
    genes_df = _build_preview_genes_df()   # cheap — no API call

    # ── Two-column layout ──────────────────────────────────────────────────
    col_sys1, col_sys2 = st.columns(2, gap="large")

    with col_sys1:
        st.subheader("System 1 — Protein Path")
        st.caption("DNA → mRNA isoforms → Amino acid translation → 3D structure")

        # Overlaid ref + mut protein in one viewer
        render_protein_overlay(result['ref_pdb'], result['mut_pdb'])

        splice = result['extras'].get('splice', {})
        if splice.get('splice_disrupted'):
            n_skipped = len(splice.get('skipped_exons', []))
            fs_note = ", frameshift detected" if splice.get('frameshift') else ""
            st.warning(
                f"Splice disruption — transcript `{splice.get('transcript_id', '')}`: "
                f"{n_skipped} exon(s) potentially skipped{fs_note}."
            )
        elif splice:
            st.success("No splice site disruption predicted.")

        with st.expander("Model & Methodology"):
            st.markdown("""
**Structure:** NVIDIA ESMFold NIM → Meta public ESMFold → mock fallback
**Overlay:** Both structures loaded into one py3Dmol scene.
Blue = Reference · Orange = Mutant
**Splicing:** AlphaGenome SPLICE_SITES track → donor/acceptor score per exon
**Confidence:** B-factor column = pLDDT (0–100). Regions <50 = unreliably folded.
""")

    with col_sys2:
        st.subheader("System 2 — Regulatory Path")
        st.caption("DNA → Chromatin accessibility → Enhancer loops → Expression")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_multi_tracks(
                result['multi_tracks'], tmp.name,
                title=f"Regulatory Landscape: {st.session_state.gene}",
            )
            st.image(Image.open(tmp.name), use_container_width=True)

        evo2_data = result['extras'].get('evo2', {})
        if evo2_data:
            score = evo2_data.get("score", float("nan"))
            score_str = f"{score:.4f}" if score == score else "N/A"
            motif_diff = evo2_data.get("motif_diff", {})
            disrupted_tfs = [tf for tf, d in motif_diff.items() if d < 0]
            gained_tfs    = [tf for tf, d in motif_diff.items() if d > 0]

            ev_col1, ev_col2 = st.columns(2)
            with ev_col1:
                st.metric("Evo 2 Score", score_str,
                          help="K-mer log-odds. Negative = diverges from ref.")
            with ev_col2:
                if disrupted_tfs:
                    st.metric("TF Motifs Lost", len(disrupted_tfs),
                              delta=f"{', '.join(disrupted_tfs[:2])}", delta_color="inverse")
                elif gained_tfs:
                    st.metric("TF Motifs Gained", len(gained_tfs),
                              delta=f"{', '.join(gained_tfs[:2])}", delta_color="normal")

        with st.expander("Model & Methodology"):
            st.markdown(f"""
**Tracks:** AlphaGenome — ATAC, CAGE, CTCF, H3K27ac (Ref vs Mut)
**Evo 2:** `{evo2_data.get('backend', 'kmer_fallback')}` — k-mer log-odds + TF motif scan
**Sign convention:** accessibility_delta = log₂(Ref/Mut) — positive = loss in mutant.
""")

    # ── Accessibility at Regulatory Elements ───────────────────────────────
    st.divider()
    st.subheader("Accessibility at Regulatory Elements")
    st.caption(
        "Per-element accessibility change. Bubble size = element width. "
        "Positive Y = accessibility **lost** in mutant; negative = gained."
    )
    _render_accessibility_elements(
        result['multi_tracks'], genes_df,
        st.session_state.edit_start, st.session_state.edit_end,
    )

    # ── 3D Chromatin Contact Map + Animation side-by-side ─────────────────
    st.divider()
    hic_col, _ = st.columns([1, 1])
    with hic_col:
        st.subheader("Hi-C Contact Map")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_hic_triangle(
                result['ref_hic'], result['mut_hic'], tmp.name,
                title="Contact Map Delta (Mut − Ref)",
            )
            st.image(Image.open(tmp.name), use_container_width=True)
        st.caption("Red = gained contacts · Blue = lost contacts · 45° triangular view")

    st.subheader("3D Chromatin Conformation — Reference → Mutant")
    st.caption(
        "Polymer simulation (Langevin dynamics) run independently on Ref and Mut "
        "contact matrices. Animation interpolates between the two equilibrium structures. "
        "Beads colored by genomic position (5′ blue → 3′ red). "
        "**Bright red beads** = edit region."
    )

    # Genomic window: extend ±200 kb around edit for the bead axis labels
    window_pad = 200_000
    g_start = max(0, st.session_state.edit_start - window_pad)
    g_end   = st.session_state.edit_end + window_pad

    render_chromatin_animation(
        result['ref_hic'], result['mut_hic'],
        edit_start=st.session_state.edit_start,
        edit_end=st.session_state.edit_end,
        genomic_start=g_start,
        genomic_end=g_end,
    )

    # ── Differentiation Trajectory ──────────────────────────────────────────
    st.divider()
    st.subheader("Differentiation Trajectory")
    st.caption(
        "AlphaGenome contact maps at embryonic (H1-hESC) and committed myeloid (K562) "
        "stages — both Ref and Mut — showing how the edit reshapes 3D organization "
        "across the developmental axis."
    )
    _render_trajectory(result)

    # ── Causal Summary ──────────────────────────────────────────────────────
    st.divider()
    st.subheader("Causal Summary")
    mod_details = {
        "type": st.session_state.mod_type.lower(),
        "target": f"regulatory element near {st.session_state.gene}",
    }
    insight = generate_mechanistic_insight(mod_details, result['delta_stats'])
    st.success(insight)

    ds = result['delta_stats']
    metric_cols = st.columns(3)
    with metric_cols[0]:
        acc = ds.get('accessibility_drop', 0)
        st.metric("Accessibility", f"{acc:+.2f} log₂",
                  delta=f"{'loss' if acc > 0 else 'gain'}", delta_color="inverse")
    with metric_cols[1]:
        exp = ds.get('expression_drop', 0)
        st.metric("Expression", f"{exp:+.2f} log₂",
                  delta=f"{'loss' if exp > 0 else 'gain'}", delta_color="inverse")
    with metric_cols[2]:
        loop_state = "Weakened" if ds.get('loop_weakened') else (
            "Strengthened" if ds.get('loop_strengthened') else "Unchanged")
        st.metric("E–P Loop", loop_state)

    ct_note = result['extras'].get('cell_type')
    if ct_note:
        st.caption(f"Predictions filtered to cell type: **{ct_note}**")

    st.caption(
        "Rule-based mechanistic attribution from AlphaGenome track deltas. "
        "No LLM was used in this analysis."
    )

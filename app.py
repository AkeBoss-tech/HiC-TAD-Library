import streamlit as st
import numpy as np
import os
import tempfile
import matplotlib.pyplot as plt
import pandas as pd
from typing import Optional
from PIL import Image
from dotenv import load_dotenv
import plotly.express as px
import plotly.graph_objects as go

load_dotenv()

from src.hg_dt.viz.tracks_plotter import plot_tracks, plot_multi_tracks
from src.hg_dt.viz.hic_plotter import plot_hic_triangle
from src.hg_dt.viz.protein_viz import render_protein_overlay
from src.hg_dt.viz.chromatin_3d import render_chromatin_animation
from src.hg_dt.viz.browser import render_genome_scroller
from src.hg_dt.analyze.attribution import generate_mechanistic_insight
from src.hg_dt.analyze.deltas import (
    accessibility_delta, expression_delta,
    find_silenced_elements, find_distal_loops,
)
from src.hg_dt.models.alphagenome import CELL_TYPES
from src.hg_dt.models.evo2 import Evo2Client

st.set_page_config(page_title="HG-DT: Human Genome Digital Twin", layout="wide")

# ---------------------------------------------------------------------------
# Gene coordinate registry (hg38) — used for quick lookup only
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
# Sidebar — API status + settings
# ---------------------------------------------------------------------------
with st.sidebar:
    st.header("Pipeline Settings")

    ag_key  = bool(os.getenv("ALPHA_GENOME_API_KEY"))
    nv_key  = bool(os.getenv("NVIDIA_API_KEY") or os.getenv("NVIDIA_ESM_FOLD_API_KEY"))

    st.markdown("**Data sources**")
    st.caption("Sequence: UCSC hg38 REST API (no key needed)")
    st.caption("Gene annotations: UCSC RefSeq + ENCODE cCREs")
    st.caption(f"AlphaGenome: {'✅ key set' if ag_key else '❌ ALPHA_GENOME_API_KEY missing'}")
    st.caption(f"ESMFold:     {'✅ NVIDIA key set' if nv_key else '⚠ using Meta public API'}")
    st.caption(f"Evo 2:       **{Evo2Client().backend}**")

    if not ag_key:
        st.error("Set ALPHA_GENOME_API_KEY in your .env file to enable predictions.")

    st.divider()
    cell_type_options = ["None (all cell types)"] + sorted(CELL_TYPES.keys())
    cell_type_sel = st.selectbox("Cell type (AlphaGenome)", cell_type_options, index=0)
    cell_type = None if cell_type_sel.startswith("None") else cell_type_sel

    st.divider()
    if st.button("Start Over"):
        for key in ["step", "chrom", "gene", "edit_start", "edit_end",
                    "mod_type", "pipeline_result", "locus_annotations"]:
            if key in st.session_state:
                del st.session_state[key]
        st.rerun()

# ---------------------------------------------------------------------------
# Session state defaults
# ---------------------------------------------------------------------------
for key, default in [
    ("step",               1),
    ("chrom",              "chr1"),
    ("gene",               "TAL1"),
    ("edit_start",         47210000),
    ("edit_end",           47215000),   # default = 5 kb deletion
    ("mod_type",           "Deletion"),
    ("pipeline_result",    None),
    ("locus_annotations",  None),
]:
    if key not in st.session_state:
        st.session_state[key] = default


# ---------------------------------------------------------------------------
# Core helpers — sequence manipulation
# ---------------------------------------------------------------------------
def _apply_edit(ref_seq: str, rel_start: int, rel_end: int,
                mod_type: str, insert_seq: str = "") -> str:
    """Apply a genomic edit to a reference sequence string."""
    t = mod_type.lower()
    if t == "deletion":
        return ref_seq[:rel_start] + ref_seq[rel_end:]
    if t == "insertion":
        payload = insert_seq if insert_seq else "A" * (rel_end - rel_start)
        return ref_seq[:rel_start] + payload + ref_seq[rel_start:]
    if t == "snp":
        pos = rel_start
        orig = ref_seq[pos]
        alt  = next(b for b in ("A", "T", "C", "G") if b != orig)
        return ref_seq[:pos] + alt + ref_seq[pos + 1:]
    return ref_seq  # "None"


# ---------------------------------------------------------------------------
# Pipeline — real data, no mocks
# ---------------------------------------------------------------------------
def _extract_1d(track_data, index: int = 0) -> np.ndarray:
    vals = np.array(track_data.values)
    return vals if vals.ndim == 1 else vals[:, index]


@st.cache_data(show_spinner="Fetching hg38 sequence from UCSC…")
def _fetch_sequences(chrom: str, edit_start: int, edit_end: int,
                     mod_type: str) -> tuple:
    """
    Fetch a 1 Mb window of hg38 from the UCSC REST API and apply the edit.
    Returns (ref_seq, mut_seq, win_start, win_end).
    """
    from src.hg_dt.data.sequence_fetcher import fetch_hg38_sequence

    center    = (edit_start + edit_end) // 2
    win_start = max(0, center - 500_000)
    win_end   = win_start + 1_000_000
    ref_seq   = fetch_hg38_sequence(chrom, win_start, win_end)

    rel_start = edit_start - win_start
    rel_end   = edit_end   - win_start
    mut_seq   = _apply_edit(ref_seq, rel_start, rel_end, mod_type)

    return ref_seq, mut_seq, win_start, win_end


@st.cache_data(show_spinner="Running AlphaGenome predictions…")
def _run_alphagenome(ref_seq: str, mut_seq: str,
                     cell_type: Optional[str]) -> tuple:
    """
    Run AlphaGenome on ref and mut sequences.
    Returns (ref_out, mut_out) AlphaGenome output objects.
    Raises RuntimeError if ALPHA_GENOME_API_KEY is not set.
    """
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector

    if not os.getenv("ALPHA_GENOME_API_KEY"):
        raise RuntimeError(
            "ALPHA_GENOME_API_KEY is not set. "
            "Add it to your .env file to run AlphaGenome predictions."
        )

    connector = AlphaGenomeConnector()
    outputs   = ["ATAC", "CAGE", "CHIP_TF", "CONTACT_MAPS", "SPLICE_SITES"]
    ref_out = connector.predict_sequence(ref_seq, organism="HUMAN",
                                          requested_outputs=outputs,
                                          cell_type=cell_type)
    mut_out = connector.predict_sequence(mut_seq, organism="HUMAN",
                                          requested_outputs=outputs,
                                          cell_type=cell_type)
    return ref_out, mut_out


@st.cache_data(show_spinner="Predicting protein structure (ESMFold)…")
def _run_protein(ref_seq: str, mut_seq: str, win_start: int,
                 edit_start: int, edit_end: int, mod_type: str,
                 chrom: str) -> tuple:
    """
    Fetch gene annotations from UCSC, extract mRNA, translate, fold with ESMFold.
    Returns (ref_pdb, mut_pdb, splice_info).
    """
    from src.hg_dt.data.sequence_fetcher import fetch_gene_annotations
    from src.hg_dt.translate.transcript import predict_isoforms
    from src.hg_dt.translate.translator import compare_translation
    from src.hg_dt.models.protein import ProteinFolder

    # Fetch real gene models from UCSC RefSeq
    genes_df = fetch_gene_annotations(chrom, win_start, win_start + 1_000_000)
    if genes_df.empty:
        return "", "", {}

    mut_shift = -(edit_end - edit_start) if mod_type.lower() == "deletion" else 0
    isoforms  = predict_isoforms(
        genes_df, ref_seq, mut_seq, win_start,
        edit_start, edit_end, mut_shift,
    )

    folder     = ProteinFolder()
    ref_pdb = mut_pdb = ""
    splice_info: dict = {}

    for tid, iso in isoforms.items():
        trans = compare_translation(iso["ref_mrna"], iso["mut_mrna"])
        if trans["ref_aa"]:
            ref_pdb = folder.predict_structure(trans["ref_aa"])["pdb"]
        if trans["mut_aa"]:
            mut_pdb = folder.predict_structure(trans["mut_aa"])["pdb"]
        splice_info = {
            "transcript_id": tid,
            "frameshift":      iso["frameshift"],
            "splice_disrupted": iso["splice_disrupted"],
            "skipped_exons":   iso["skipped_exons"],
        }
        break   # first transcript only

    return ref_pdb, mut_pdb, splice_info


@st.cache_data(show_spinner="Running Evo 2 variant scoring…")
def _run_evo2(ref_seq: str, mut_seq: str) -> tuple:
    window  = min(2048, len(ref_seq))
    center  = len(ref_seq) // 2
    ref_win = ref_seq[center - window//2 : center + window//2]
    mut_win = mut_seq[center - window//2 : center + window//2]
    evo2    = Evo2Client()
    return evo2.score_variant(ref_win, mut_win), evo2.scan_sequence(ref_win, window=256)


@st.cache_data(show_spinner="Running differentiation trajectory…")
def _run_trajectory(ref_seq: str, mut_seq: str) -> dict:
    """
    Run AlphaGenome at H1-hESC and K562 for both sequences.
    Returns 4 contact maps.
    """
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector

    if not os.getenv("ALPHA_GENOME_API_KEY"):
        raise RuntimeError("ALPHA_GENOME_API_KEY not set.")

    connector = AlphaGenomeConnector()
    maps = {}
    for label, seq, ct in [
        ("h1_ref",   ref_seq, "H1-hESC"),
        ("h1_mut",   mut_seq, "H1-hESC"),
        ("k562_ref", ref_seq, "K562"),
        ("k562_mut", mut_seq, "K562"),
    ]:
        out = connector.predict_sequence(seq, requested_outputs=["CONTACT_MAPS"],
                                          cell_type=ct)
        maps[label] = np.array(out.contact_maps.values[:, :, 0])
    return maps


def _run_pipeline() -> dict:
    """
    Orchestrate the full pipeline and return a result dict.
    Raises on any error — no mock fallback.
    """
    chrom      = st.session_state.chrom
    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end
    mod_type   = st.session_state.mod_type

    # 1. Real sequences from UCSC
    ref_seq, mut_seq, win_start, win_end = _fetch_sequences(
        chrom, edit_start, edit_end, mod_type
    )

    # 2. AlphaGenome predictions
    ref_out, mut_out = _run_alphagenome(ref_seq, mut_seq, cell_type)

    # 3. Extract tracks
    ref_atac = _extract_1d(ref_out.atac)
    mut_atac = _extract_1d(mut_out.atac)
    ref_cage = _extract_1d(ref_out.cage)
    mut_cage = _extract_1d(mut_out.cage)
    ref_ctcf = _extract_1d(ref_out.chip_tf)
    mut_ctcf = _extract_1d(mut_out.chip_tf)
    ref_hic  = np.array(ref_out.contact_maps.values[:, :, 0])
    mut_hic  = np.array(mut_out.contact_maps.values[:, :, 0])

    multi_tracks = {
        "ATAC": (ref_atac, mut_atac),
        "CAGE": (ref_cage, mut_cage),
        "CTCF": (ref_ctcf, mut_ctcf),
    }

    # 4. Delta statistics
    a_delta        = accessibility_delta(ref_atac, mut_atac)
    e_delta        = expression_delta(ref_cage, mut_cage)
    silenced       = find_silenced_elements(ref_atac, mut_atac)
    loop_weakened  = len(find_distal_loops(mut_hic)) < len(find_distal_loops(ref_hic)) * 0.8
    loop_gained    = len(find_distal_loops(mut_hic)) > len(find_distal_loops(ref_hic)) * 1.2
    acc_drop       = float(np.mean(a_delta[silenced])) if silenced else float(np.mean(a_delta[a_delta > 0]) or 0)
    exp_drop       = float(np.mean(e_delta[e_delta > 0]) or 0)

    delta_stats = {
        "loop_weakened":    loop_weakened,
        "loop_strengthened": loop_gained,
        "accessibility_drop": acc_drop,
        "expression_drop":    exp_drop,
    }

    # 5. Protein (parallel to tracks — can fail gracefully)
    try:
        ref_pdb, mut_pdb, splice_info = _run_protein(
            ref_seq, mut_seq, win_start, edit_start, edit_end, mod_type, chrom
        )
    except Exception as exc:
        ref_pdb = mut_pdb = ""
        splice_info = {"error": str(exc)}

    # 6. Evo 2
    evo2_result, evo2_scan = _run_evo2(ref_seq, mut_seq)

    return {
        "ref_seq":      ref_seq,
        "mut_seq":      mut_seq,
        "win_start":    win_start,
        "ref_atac":     ref_atac,
        "mut_atac":     mut_atac,
        "ref_hic":      ref_hic,
        "mut_hic":      mut_hic,
        "ref_pdb":      ref_pdb,
        "mut_pdb":      mut_pdb,
        "multi_tracks": multi_tracks,
        "delta_stats":  delta_stats,
        "extras": {
            "evo2":   evo2_result,
            "evo2_scan": evo2_scan,
            "splice": splice_info,
            "cell_type": cell_type,
        },
    }


# ---------------------------------------------------------------------------
# Helper: fetch real locus annotations for Step 2
# ---------------------------------------------------------------------------
@st.cache_data(show_spinner="Fetching gene annotations from UCSC…")
def _fetch_locus_annotations(chrom: str, edit_start: int,
                              edit_end: int) -> pd.DataFrame:
    from src.hg_dt.data.sequence_fetcher import build_locus_annotations
    return build_locus_annotations(chrom, edit_start, edit_end, window=500_000)


# ---------------------------------------------------------------------------
# Helper: render per-element accessibility bubble chart
# ---------------------------------------------------------------------------
def _render_accessibility_elements(multi_tracks: dict, genes_df: pd.DataFrame,
                                    edit_start: int, edit_end: int):
    atac_ref, atac_mut = multi_tracks.get("ATAC", (None, None))
    if atac_ref is None or genes_df is None or genes_df.empty:
        return

    n           = len(atac_ref)
    window_span = max(edit_end - edit_start, 1)
    rows = []

    for _, row in genes_df.iterrows():
        mid  = (row["Start"] + row["End"]) / 2
        size = max(row["End"] - row["Start"], 1)
        idx  = int((mid - edit_start) / window_span * n)
        idx  = max(0, min(n - 1, idx))
        i0, i1 = max(0, idx - 2), min(n, idx + 3)

        ref_sig = float(np.mean(atac_ref[i0:i1]))
        mut_sig = float(np.mean(atac_mut[i0:i1]))
        delta   = float(np.log2((ref_sig + 1e-6) / (mut_sig + 1e-6)))

        rows.append({
            "Element":           row["Name"],
            "Type":              row["Type"],
            "Position":          int(mid),
            "Element size (bp)": int(size),
            "Ref ATAC":          round(ref_sig, 3),
            "Mut ATAC":          round(mut_sig, 3),
            "Δ log₂(Ref/Mut)":   round(delta, 3),
        })

    if not rows:
        return

    df  = pd.DataFrame(rows)
    fig = px.scatter(
        df, x="Position", y="Δ log₂(Ref/Mut)",
        size="Element size (bp)", color="Type", text="Element",
        hover_data=["Ref ATAC", "Mut ATAC"],
        color_discrete_map={"Gene": "#4a9eff", "Enhancer": "#ff9a3c",
                             "CTCF": "#2ecc71", "Promoter": "#a855f7",
                             "cCRE": "#f43f5e", "Insulator": "#14b8a6"},
        size_max=45,
    )
    fig.add_hline(y=0, line_dash="dash", line_color="rgba(180,180,180,0.5)")
    fig.add_vrect(x0=edit_start, x1=edit_end,
                  fillcolor="rgba(255,80,80,0.08)",
                  line_color="rgba(255,80,80,0.4)",
                  annotation_text="Edit", annotation_position="top left")
    fig.update_traces(textposition="top center", textfont_size=10)
    fig.update_layout(
        height=320,
        xaxis=dict(tickformat=",.0s", title="Genomic position (hg38)"),
        yaxis=dict(title="Accessibility Δ log₂ (Ref / Mut)", zeroline=False),
        plot_bgcolor="white", paper_bgcolor="white",
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        margin=dict(t=30, b=40, l=60, r=20),
    )
    st.plotly_chart(fig, use_container_width=True)


# ---------------------------------------------------------------------------
# Helper: render differentiation trajectory
# ---------------------------------------------------------------------------
def _render_trajectory(result: dict):
    try:
        maps = _run_trajectory(result["ref_seq"], result["mut_seq"])
    except Exception as exc:
        st.error(f"Differentiation trajectory failed: {exc}")
        return

    titles = ["H1-hESC  Ref", "H1-hESC  Mut", "K562  Ref", "K562  Mut"]
    keys   = ["h1_ref", "h1_mut", "k562_ref", "k562_mut"]
    vmax   = max(np.nanmax(np.abs(maps[k])) for k in keys if maps[k].size > 0)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for ax, key, title in zip(axes, keys, titles):
        im = ax.imshow(maps[key], cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="Contact frequency")
    fig.text(0.5, -0.02,
             "← Embryonic stem cell (H1-hESC)   ·   Committed myeloid (K562) →",
             ha="center", fontsize=10, color="#444")
    plt.tight_layout()
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        fig.savefig(tmp.name, dpi=150, bbox_inches="tight")
        plt.close(fig)
        st.image(Image.open(tmp.name), use_container_width=True)

    with st.expander("Model & Methodology"):
        st.markdown("""
**Model:** AlphaGenome (DeepMind) with cell-type-specific ontology terms
**Cell types:** H1-hESC (`CLO:0037111`) = embryonic pluripotent stem cell; K562 (`CLO:0007050`) = committed CML myeloid line
**Interpretation:** Contact maps predicted from real hg38 sequence for each cell type. Changes between Ref and Mut show how the edit disrupts 3D genome organization across the developmental axis.
""")


# ---------------------------------------------------------------------------
# Step progress bar
# ---------------------------------------------------------------------------
st.title("HG-DT: Human Genome Digital Twin")

step_labels = ["1 · Find Locus", "2 · What's Here?", "3 · Results"]
prog_cols   = st.columns(3)
for i, (col, label) in enumerate(zip(prog_cols, step_labels)):
    active = (i + 1 == st.session_state.step)
    col.markdown(
        f"**:blue[{label}]**" if active
        else f"<span style='color:#aaa'>{label}</span>",
        unsafe_allow_html=True,
    )

st.divider()


# ===========================================================================
# STEP 1 — Find Your Locus
# ===========================================================================
if st.session_state.step == 1:
    st.header("Step 1 — Find Your Locus")
    st.caption("Search for a gene or enter coordinates, then define your edit.")

    col_search, col_jump = st.columns([3, 1])
    with col_search:
        search_input = st.text_input(
            "Gene symbol or chr:start-end",
            value=st.session_state.gene,
            placeholder="e.g. TAL1  or  chr1:47210000-47260000",
        )
    with col_jump:
        st.write("")
        if st.button("Look up", use_container_width=True):
            g = search_input.upper().strip()
            if g in GENE_COORDS:
                c, s, e = GENE_COORDS[g]
                st.session_state.chrom      = c
                st.session_state.edit_start = s
                st.session_state.edit_end   = e
                st.session_state.gene       = g
                st.rerun()
            elif ":" in search_input and "-" in search_input:
                try:
                    chrom_p, range_p = search_input.replace(",", "").split(":")
                    s, e = range_p.split("-")
                    st.session_state.chrom      = chrom_p.strip()
                    st.session_state.edit_start = int(s.strip())
                    st.session_state.edit_end   = int(e.strip())
                    st.session_state.gene       = chrom_p.strip()
                    st.rerun()
                except Exception:
                    st.error("Could not parse. Use chr1:47210000-47260000")
            else:
                st.warning(f"'{search_input}' not in registry — enter coordinates manually.")

    st.write("")
    col_mod, col_chrom = st.columns(2)
    with col_mod:
        mod_options = ["Deletion", "Insertion", "SNP", "None"]
        mod_type = st.selectbox(
            "Modification type",
            mod_options,
            index=mod_options.index(st.session_state.mod_type)
                  if st.session_state.mod_type in mod_options else 0,
        )
    with col_chrom:
        chrom = st.text_input("Chromosome", value=st.session_state.chrom)

    col_s, col_e = st.columns(2)
    with col_s:
        edit_start = st.number_input("Edit start (bp)", value=st.session_state.edit_start,
                                     step=1000)
    with col_e:
        edit_end = st.number_input("Edit end (bp)", value=st.session_state.edit_end,
                                   step=1000)

    edit_size = int(edit_end) - int(edit_start)
    if edit_size > 0:
        st.caption(
            f"Edit: **{edit_size:,} bp {mod_type.lower()}** "
            f"({edit_size / 1000:.1f} kb)"
        )
        # Quick-select common deletion sizes
        if mod_type == "Deletion":
            st.caption("Quick-select deletion size:")
            qcols = st.columns(5)
            sizes = [1_000, 5_000, 10_000, 50_000, 100_000]
            labels = ["1 kb", "5 kb", "10 kb", "50 kb", "100 kb"]
            for qcol, size, lbl in zip(qcols, sizes, labels):
                if qcol.button(lbl, use_container_width=True):
                    st.session_state.edit_end = int(edit_start) + size
                    st.rerun()
    elif edit_size <= 0:
        st.error("Edit end must be greater than edit start.")

    st.write("")
    if st.button("Preview Locus →", type="primary", disabled=(edit_size <= 0)):
        st.session_state.chrom      = chrom
        st.session_state.edit_start = int(edit_start)
        st.session_state.edit_end   = int(edit_end)
        st.session_state.mod_type   = mod_type
        st.session_state.locus_annotations = None   # force refresh
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

    # Fetch real annotations (cached)
    if st.session_state.locus_annotations is None:
        with st.spinner("Fetching gene annotations from UCSC…"):
            try:
                annot = _fetch_locus_annotations(
                    st.session_state.chrom, edit_start, edit_end
                )
                st.session_state.locus_annotations = annot
            except Exception as exc:
                st.warning(f"UCSC annotation fetch failed: {exc}")
                st.session_state.locus_annotations = pd.DataFrame(
                    columns=["Name", "Type", "Start", "End"]
                )

    genes_df = st.session_state.locus_annotations

    render_genome_scroller(
        genes_df,
        current_pos=mid_pos,
        window_size=500_000,
        deletion_region=(edit_start, edit_end),
    )

    # Overlap table
    if not genes_df.empty:
        overlapping = genes_df[
            (genes_df["Start"] <= edit_end) & (genes_df["End"] >= edit_start)
        ][["Name", "Type", "Start", "End"]].copy()

        if not overlapping.empty:
            st.warning(f"Planned edit overlaps **{len(overlapping)}** annotated element(s):")
            overlapping["Start"] = overlapping["Start"].map("{:,}".format)
            overlapping["End"]   = overlapping["End"].map("{:,}".format)
            st.dataframe(overlapping, use_container_width=True, hide_index=True)
        else:
            st.info("Planned edit does not overlap any annotated elements in this view.")
    else:
        st.info("No annotations found (UCSC may be unavailable). You can still predict.")

    st.write("")
    if st.button("Predict Effect →", type="primary"):
        st.session_state.step            = 3
        st.session_state.pipeline_result = None
        st.rerun()


# ===========================================================================
# STEP 3 — Results
# ===========================================================================
elif st.session_state.step == 3:
    if st.button("← Back to Preview"):
        st.session_state.step = 2
        st.rerun()

    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end

    st.header("Step 3 — Predicted Effect")
    st.caption(
        f"**{st.session_state.gene}** · "
        f"{st.session_state.chrom}:{edit_start:,}–{edit_end:,} · "
        f"{edit_end - edit_start:,} bp {st.session_state.mod_type.lower()}"
        + (f" · Cell type: {cell_type}" if cell_type else "")
    )

    # Run pipeline (cached in session_state)
    if st.session_state.pipeline_result is None:
        try:
            st.session_state.pipeline_result = _run_pipeline()
        except Exception as exc:
            st.error(f"Pipeline error: {exc}")
            st.stop()

    result   = st.session_state.pipeline_result
    genes_df = st.session_state.locus_annotations or pd.DataFrame(
        columns=["Name", "Type", "Start", "End"]
    )

    # ── Two-column layout ──────────────────────────────────────────────────
    col_sys1, col_sys2 = st.columns(2, gap="large")

    with col_sys1:
        st.subheader("System 1 — Protein Path")
        st.caption("DNA → mRNA isoforms → Amino acid translation → 3D structure")

        render_protein_overlay(result["ref_pdb"], result["mut_pdb"])

        splice = result["extras"].get("splice", {})
        if splice.get("splice_disrupted"):
            n_skipped = len(splice.get("skipped_exons", []))
            st.warning(
                f"Splice disruption: transcript `{splice.get('transcript_id', '')}` — "
                f"{n_skipped} exon(s) skipped"
                + (", frameshift" if splice.get("frameshift") else "") + "."
            )
        elif splice and not splice.get("error"):
            st.success("No splice site disruption predicted.")
        elif splice.get("error"):
            st.info(f"Protein prediction note: {splice['error']}")

        with st.expander("Model & Methodology"):
            st.markdown("""
**Structure:** NVIDIA ESMFold NIM → Meta public ESMFold (free)
**Overlay:** Blue = Reference · Orange = Mutant (both in one py3Dmol scene)
**Splicing:** AlphaGenome SPLICE_SITES track → donor/acceptor score per exon
**Gene models:** UCSC hg38 RefSeq (live REST API)
**Confidence:** B-factor column = pLDDT (0–100). Regions <50 = unreliably folded.
""")

    with col_sys2:
        st.subheader("System 2 — Regulatory Path")
        st.caption("DNA → Chromatin accessibility → Enhancer loops → Expression")

        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_multi_tracks(
                result["multi_tracks"], tmp.name,
                title=f"AlphaGenome tracks: {st.session_state.gene}",
            )
            st.image(Image.open(tmp.name), use_container_width=True)

        evo2_data = result["extras"].get("evo2", {})
        if evo2_data:
            score      = evo2_data.get("score", float("nan"))
            motif_diff = evo2_data.get("motif_diff", {})
            disrupted  = [tf for tf, d in motif_diff.items() if d < 0]
            gained     = [tf for tf, d in motif_diff.items() if d > 0]
            ev1, ev2   = st.columns(2)
            with ev1:
                st.metric("Evo 2 Score", f"{score:.4f}" if score == score else "N/A",
                          help="K-mer log-odds. Negative = alt diverges from ref.")
            with ev2:
                if disrupted:
                    st.metric("TF Motifs Lost", len(disrupted),
                              delta=", ".join(disrupted[:2]), delta_color="inverse")
                elif gained:
                    st.metric("TF Motifs Gained", len(gained),
                              delta=", ".join(gained[:2]), delta_color="normal")

        with st.expander("Model & Methodology"):
            st.markdown(f"""
**Tracks:** AlphaGenome — ATAC, CAGE, CTCF (Ref vs Mut); sequence from UCSC hg38
**Evo 2:** `{evo2_data.get('backend', 'kmer_fallback')}` — k-mer log-odds + CTCF/SP1/GATA1 motif scan
**Sign convention:** Δ = log₂(Ref/Mut) — positive = loss in mutant
""")

    # ── Accessibility at Regulatory Elements ───────────────────────────────
    st.divider()
    st.subheader("Accessibility at Regulatory Elements")
    st.caption(
        "Per-element accessibility change from real ENCODE cCREs + RefSeq genes (UCSC). "
        "Bubble size = element width. Positive Y = accessibility **lost** in mutant."
    )
    _render_accessibility_elements(
        result["multi_tracks"], genes_df, edit_start, edit_end
    )

    # ── Hi-C + 3D animation ────────────────────────────────────────────────
    st.divider()
    hic_col, _ = st.columns([1, 1])
    with hic_col:
        st.subheader("Hi-C Contact Map")
        with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
            plot_hic_triangle(
                result["ref_hic"], result["mut_hic"], tmp.name,
                title="Contact Map Delta (Mut − Ref)",
            )
            st.image(Image.open(tmp.name), use_container_width=True)
        st.caption("Red = gained contacts · Blue = lost contacts · triangular view")

    st.subheader("3D Chromatin Conformation — Reference → Mutant")
    st.caption(
        "Polymer simulation (Langevin dynamics) on AlphaGenome contact maps. "
        "Press ▶ to watch chromatin reorganize. "
        "Rainbow = genomic position (5′→3′). Red beads = edit region."
    )
    window_pad = 200_000
    render_chromatin_animation(
        result["ref_hic"], result["mut_hic"],
        edit_start=edit_start,
        edit_end=edit_end,
        genomic_start=max(0, edit_start - window_pad),
        genomic_end=edit_end + window_pad,
    )

    # ── Differentiation Trajectory ─────────────────────────────────────────
    st.divider()
    st.subheader("Differentiation Trajectory")
    st.caption(
        "AlphaGenome contact maps at H1-hESC (embryonic stem) → K562 (committed myeloid). "
        "Shows how the edit reshapes 3D organization across the developmental axis."
    )
    _render_trajectory(result)

    # ── Causal Summary ──────────────────────────────────────────────────────
    st.divider()
    st.subheader("Causal Summary")
    mod_details = {
        "type":   st.session_state.mod_type.lower(),
        "target": f"regulatory element near {st.session_state.gene}",
    }
    insight = generate_mechanistic_insight(mod_details, result["delta_stats"])
    st.success(insight)

    ds        = result["delta_stats"]
    mc1, mc2, mc3 = st.columns(3)
    with mc1:
        acc = ds.get("accessibility_drop", 0)
        st.metric("Accessibility", f"{acc:+.2f} log₂",
                  delta="loss" if acc > 0 else "gain", delta_color="inverse")
    with mc2:
        exp = ds.get("expression_drop", 0)
        st.metric("Expression", f"{exp:+.2f} log₂",
                  delta="loss" if exp > 0 else "gain", delta_color="inverse")
    with mc3:
        loop_state = ("Weakened"    if ds.get("loop_weakened")    else
                      "Strengthened" if ds.get("loop_strengthened") else "Unchanged")
        st.metric("E–P Loop", loop_state)

    ct_note = result["extras"].get("cell_type")
    if ct_note:
        st.caption(f"Predictions filtered to cell type: **{ct_note}**")
    st.caption(
        "Mechanistic attribution from AlphaGenome track deltas (log₂ Ref/Mut). "
        "Sequence from UCSC hg38. No hardcoded or mock data."
    )

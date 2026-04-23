import streamlit as st
import numpy as np
import os
import tempfile
import json
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
from src.hg_dt.viz.chromatin_3d import render_chromatin_animation, run_chromatin_sim
from src.hg_dt.viz.browser import render_genome_scroller
from src.hg_dt.analyze.attribution import generate_mechanistic_insight
from src.hg_dt.analyze.deltas import (
    accessibility_delta, expression_delta,
    find_silenced_elements, find_distal_loops,
)
from src.hg_dt.models.alphagenome import CELL_TYPES
from src.hg_dt.models.evo2 import Evo2Client
from src.hg_dt.data.reference_paths import (
    GENCODE_GTF,
    CCRE_BED,
    hg38_fasta_available,
    gencode_gtf_available,
    ccre_bed_available,
    preferred_hg38_fasta,
)

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
    seq_src = (
        "local hg38 FASTA (`data/references/hg38.fa`)"
        if hg38_fasta_available()
        else "UCSC hg38 REST API (no key needed)"
    )
    st.caption(f"Sequence: {seq_src}")
    st.caption("Gene annotations: UCSC RefSeq + ENCODE cCREs (API)")
    if gencode_gtf_available():
        st.caption(f"Local GENCODE GTF: ✅ `{GENCODE_GTF}`")
    if ccre_bed_available():
        st.caption(f"Local SCREEN cCRE BED: ✅ `{CCRE_BED}`")
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
                    "mod_type", "insert_seq", "snp_alt",
                    "pipeline_result", "locus_annotations"]:
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
    ("insert_seq",         ""),          # insertion: user bases (length = edit span)
    ("snp_alt",            "A"),         # SNP alternate base (validated vs ref at run time)
    ("pipeline_result",    None),
    ("locus_annotations",  None),
]:
    if key not in st.session_state:
        st.session_state[key] = default


# ---------------------------------------------------------------------------
# Pipeline — real data, no mocks
# ---------------------------------------------------------------------------
def _extract_1d(track_data, index: int = 0) -> np.ndarray:
    vals = np.array(track_data.values)
    return vals if vals.ndim == 1 else vals[:, index]


@st.cache_data(show_spinner="Fetching hg38 sequence…")
def _fetch_sequences(
    chrom: str,
    edit_start: int,
    edit_end: int,
    mod_type: str,
    insert_seq: str = "",
    snp_alt: Optional[str] = None,
    fasta_path: Optional[str] = None,
) -> tuple:
    """
    Fetch hg38 windows and apply the edit so that BOTH ref and mut are
    exactly 1,048,576 bp — required by AlphaGenome.

    See :func:`src.hg_dt.data.sequence_fetcher.fetch_ref_mut_sequences`.
    """
    from src.hg_dt.data.sequence_fetcher import fetch_ref_mut_sequences

    fp = fasta_path if fasta_path is not None else preferred_hg38_fasta()
    return fetch_ref_mut_sequences(
        chrom,
        edit_start,
        edit_end,
        mod_type,
        insert_seq=insert_seq,
        snp_alt=snp_alt,
        fasta_path=fp,
    )


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
def _run_protein(
    ref_seq: str,
    mut_seq: str,
    win_start: int,
    edit_start: int,
    edit_end: int,
    mod_type: str,
    chrom: str,
    ref_out=None,
    mut_out=None,
) -> tuple:
    """
    UCSC RefSeq models → :func:`predict_isoforms` + longest-ORF translation
    (``compare_translation``), optional SPLICE_SITES from AlphaGenome, then
    ESMFold. Falls back to CDS translation if ORF/isoform parsing is empty.
    Returns (ref_pdb, mut_pdb, splice_info).
    """
    from src.hg_dt.data.sequence_fetcher import fetch_gene_models, compute_mut_shift
    from src.hg_dt.translate.transcript import ucsc_refseq_models_to_genes_df, predict_isoforms
    from src.hg_dt.translate.translator import compare_translation
    from src.hg_dt.models.protein import ProteinFolder

    fold_cap = ProteinFolder.ESMFOLD_MAX_INPUT_AA

    def _extract_splice_1d(out):
        if out is None:
            return None
        ss = getattr(out, "splice_sites", None)
        if ss is None:
            return None
        try:
            vals = np.array(ss.values)
            return vals if vals.ndim == 1 else vals[:, 0]
        except Exception:
            return None

    win_end = win_start + len(ref_seq)
    models = fetch_gene_models(chrom, win_start, win_end)
    if not models:
        return "", "", {"error": "No gene models in window"}

    def _score_model(m):
        if m["cdsStart"] >= m["cdsEnd"]:
            return float("inf")
        if m["cdsStart"] < win_start or m["cdsEnd"] > win_end:
            return float("inf")
        gene_mid = (m["txStart"] + m["txEnd"]) / 2
        edit_mid = (edit_start + edit_end) / 2
        overlap = max(0, min(m["txEnd"], edit_end) - max(m["txStart"], edit_start))
        if overlap > 0:
            return 0
        return abs(gene_mid - edit_mid)

    scored = sorted(models, key=_score_model)
    chosen = next((m for m in scored if _score_model(m) < float("inf")), None)
    if chosen is None:
        return "", "", {"error": "No coding gene within window"}

    def _assemble_cds(seq: str, win_s: int, exon_starts, exon_ends,
                      cds_start: int, cds_end: int, strand: str) -> str:
        cds_dna = ""
        for es, ee in zip(exon_starts, exon_ends):
            s = max(es, cds_start)
            e = min(ee, cds_end)
            if s >= e:
                continue
            rel_s = s - win_s
            rel_e = e - win_s
            if rel_s < 0 or rel_e > len(seq):
                continue
            cds_dna += seq[rel_s:rel_e]
        if strand == "-":
            comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
            cds_dna = cds_dna.translate(comp)[::-1]
        return cds_dna.upper()

    def _translate_cds(dna: str) -> str:
        codon_table = {
            "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
        }
        protein = ""
        for i in range(0, len(dna) - 2, 3):
            codon = dna[i:i+3]
            aa = codon_table.get(codon, "X")
            if aa == "*":
                break
            protein += aa
        return protein

    genes_df = ucsc_refseq_models_to_genes_df(models)
    ins_len = edit_end - edit_start
    mut_shift = compute_mut_shift(mod_type, edit_start, edit_end, ins_len)
    ref_sp = _extract_splice_1d(ref_out)
    mut_sp = _extract_splice_1d(mut_out)

    isoforms = {}
    if not genes_df.empty:
        isoforms = predict_isoforms(
            genes_df, ref_seq, mut_seq, win_start,
            edit_start, edit_end, mut_shift,
            ref_splice_track=ref_sp,
            mut_splice_track=mut_sp,
        )

    tid = chosen["name"]
    iso = isoforms.get(tid) if isoforms else None
    ref_aa, mut_aa = "", ""
    comp = {}
    skipped: list = []
    splice_disrupted = False
    frameshift = False

    if iso:
        comp = compare_translation(iso["ref_mrna"], iso["mut_mrna"])
        ref_aa = comp.get("ref_aa", "") or ""
        mut_aa = comp.get("mut_aa", "") or ""
        skipped = list(iso.get("skipped_exons", []))
        splice_disrupted = bool(iso.get("splice_disrupted"))
        frameshift = bool(comp.get("is_frameshift") or iso.get("frameshift"))

    if not ref_aa:
        ref_cds = _assemble_cds(
            ref_seq, win_start,
            chosen["exon_starts"], chosen["exon_ends"],
            chosen["cdsStart"], chosen["cdsEnd"], chosen["strand"],
        )
        mut_cds = _assemble_cds(
            mut_seq, win_start,
            chosen["exon_starts"], chosen["exon_ends"],
            chosen["cdsStart"], chosen["cdsEnd"], chosen["strand"],
        )
        ref_aa = _translate_cds(ref_cds)
        mut_aa = _translate_cds(mut_cds)
        len_diff = len(mut_aa) - len(ref_aa)
        frameshift = (len_diff % 3 != 0) if mut_aa else False
        splice_disrupted = len_diff != 0 and abs(len_diff) >= 5
        skipped = []

    if not ref_aa:
        return "", "", {"error": f"Empty translation for {chosen['name2']}"}

    ref_aa_fold = ref_aa[:fold_cap]
    mut_aa_fold = mut_aa[:fold_cap] if mut_aa else ref_aa_fold

    folder = ProteinFolder()
    ref_pdb = folder.predict_structure(ref_aa_fold)["pdb"]
    mut_pdb = folder.predict_structure(mut_aa_fold)["pdb"]

    splice_info = {
        "transcript_id": tid,
        "gene": chosen["name2"],
        "ref_aa_len": len(ref_aa),
        "mut_aa_len": len(mut_aa),
        "frameshift": frameshift,
        "splice_disrupted": splice_disrupted,
        "skipped_exons": skipped,
        "truncated_to_fold": len(ref_aa) > fold_cap or len(mut_aa) > fold_cap,
        "fold_cap_aa": fold_cap,
    }
    return ref_pdb, mut_pdb, splice_info


@st.cache_data(show_spinner="Running Evo 2 variant scoring…")
def _run_evo2(ref_seq: str, mut_seq: str) -> tuple:
    window  = min(2048, len(ref_seq))
    center  = len(ref_seq) // 2
    ref_win = ref_seq[center - window//2 : center + window//2]
    mut_win = mut_seq[center - window//2 : center + window//2]
    evo2    = Evo2Client()
    return evo2.score_variant(ref_win, mut_win), evo2.scan_sequence(ref_win, window=256)


@st.cache_data(show_spinner=False)
def _run_trajectory(ref_seq: str, mut_seq: str) -> dict:
    """
    Run AlphaGenome contact maps at H1-hESC and GM12878 for ref + mut.
    Returns 4 contact maps.
    """
    import time
    from src.hg_dt.models.alphagenome import AlphaGenomeConnector

    if not os.getenv("ALPHA_GENOME_API_KEY"):
        raise RuntimeError("ALPHA_GENOME_API_KEY not set.")

    connector = AlphaGenomeConnector()
    maps = {}
    # Only H1 and GM12878 have AlphaGenome contact map tracks (K562 has none)
    labels = [
        ("h1_ref",      ref_seq, "H1",      "Ref @ H1-hESC (stem)"),
        ("h1_mut",      mut_seq, "H1",      "Mut @ H1-hESC (stem)"),
        ("gm12878_ref", ref_seq, "GM12878", "Ref @ GM12878 (committed)"),
        ("gm12878_mut", mut_seq, "GM12878", "Mut @ GM12878 (committed)"),
    ]
    for i, (label, seq, ct, display) in enumerate(labels, 1):
        st.write(f"  [{i}/4] AlphaGenome → `{display}`…")
        t = time.time()
        out = connector.predict_sequence(seq, requested_outputs=["CONTACT_MAPS"],
                                          cell_type=ct)
        cm = np.array(out.contact_maps.values)
        if cm.ndim == 3 and cm.shape[2] > 0:
            maps[label] = cm[:, :, 0]
        elif cm.ndim == 2:
            maps[label] = cm
        else:
            raise RuntimeError(
                f"No contact map tracks available for cell type {ct!r}. "
                f"AlphaGenome only supports Hi-C for: H1, GM12878, HFFc6, "
                f"HCT116, IMR-90, HeLa, HepG2."
            )
        st.write(f"    ✅ contact map {maps[label].shape}  · _{time.time() - t:.1f}s_")
    return maps


def _run_pipeline() -> dict:
    """
    Orchestrate the full pipeline and return a result dict.
    Logs each step via st.write() — call inside st.status() for a live log panel.
    Raises on any error — no mock fallback.
    """
    import time
    t0 = time.time()

    chrom      = st.session_state.chrom
    edit_start = st.session_state.edit_start
    edit_end   = st.session_state.edit_end
    mod_type   = st.session_state.mod_type
    edit_size  = edit_end - edit_start

    # ── Step 1: Fetch sequences ────────────────────────────────────────────
    st.write(
        f"**[1/5]** Fetching hg38 sequence from UCSC  "
        f"(`{chrom}:{edit_start:,}–{edit_end:,}`, {edit_size:,} bp {mod_type.lower()})…"
    )
    t1 = time.time()
    ins_raw = st.session_state.get("insert_seq", "") or ""
    snp_a = st.session_state.get("snp_alt")
    ref_seq, mut_seq, win_start, win_end = _fetch_sequences(
        chrom, edit_start, edit_end, mod_type,
        insert_seq=ins_raw,
        snp_alt=snp_a if mod_type == "SNP" else None,
    )
    st.write(
        f"  ✅ Window: `{chrom}:{win_start:,}–{win_end:,}` ({win_end - win_start:,} bp)  "
        f"· ref = {len(ref_seq):,} bp  · mut = {len(mut_seq):,} bp  "
        f"· _{time.time() - t1:.1f}s_"
    )

    # ── Step 2: AlphaGenome ────────────────────────────────────────────────
    st.write(
        f"**[2/5]** Running AlphaGenome on ref + mut sequences "
        f"(cell type: {cell_type or 'all'})…"
    )
    t2 = time.time()
    ref_out, mut_out = _run_alphagenome(ref_seq, mut_seq, cell_type)
    st.write(
        f"  ✅ Predictions complete — ATAC · CAGE · CTCF · Contact maps · Splice sites  "
        f"· _{time.time() - t2:.1f}s_"
    )

    # ── Step 3: Extract tracks & delta stats ──────────────────────────────
    st.write("**[3/5]** Extracting tracks and computing delta statistics…")
    t3 = time.time()
    ref_atac = _extract_1d(ref_out.atac)
    mut_atac = _extract_1d(mut_out.atac)
    ref_cage = _extract_1d(ref_out.cage)
    mut_cage = _extract_1d(mut_out.cage)
    ref_ctcf = _extract_1d(ref_out.chip_tf)
    mut_ctcf = _extract_1d(mut_out.chip_tf)
    def _extract_hic(out) -> np.ndarray:
        cm = np.array(out.contact_maps.values)
        if cm.ndim == 3 and cm.shape[2] > 0:
            return cm[:, :, 0]
        if cm.ndim == 2:
            return cm
        # No contact maps for this cell type — return zeros
        n = 512
        return np.zeros((n, n))

    ref_hic  = _extract_hic(ref_out)
    mut_hic  = _extract_hic(mut_out)

    multi_tracks = {
        "ATAC": (ref_atac, mut_atac),
        "CAGE": (ref_cage, mut_cage),
        "CTCF": (ref_ctcf, mut_ctcf),
    }

    a_delta       = accessibility_delta(ref_atac, mut_atac)
    e_delta       = expression_delta(ref_cage, mut_cage)
    silenced      = find_silenced_elements(ref_atac, mut_atac)
    ref_loops     = find_distal_loops(ref_hic)
    mut_loops     = find_distal_loops(mut_hic)
    loop_weakened = len(mut_loops) < len(ref_loops) * 0.8
    loop_gained   = len(mut_loops) > len(ref_loops) * 1.2
    acc_drop      = float(np.mean(a_delta[silenced])) if silenced else float(np.mean(a_delta[a_delta > 0]) or 0)
    exp_drop      = float(np.mean(e_delta[e_delta > 0]) or 0)

    delta_stats = {
        "loop_weakened":     loop_weakened,
        "loop_strengthened": loop_gained,
        "accessibility_drop": acc_drop,
        "expression_drop":    exp_drop,
    }
    st.write(
        f"  ✅ Hi-C: {ref_hic.shape[0]}×{ref_hic.shape[1]}  "
        f"· ATAC bins: {len(ref_atac):,}  "
        f"· Silenced elements: {len(silenced)}  "
        f"· Accessibility Δ = {acc_drop:+.3f}  "
        f"· Expression Δ = {exp_drop:+.3f}  "
        f"· Loops ref/mut = {len(ref_loops)}/{len(mut_loops)}  "
        f"· _{time.time() - t3:.1f}s_"
    )

    # ── Step 4: Protein folding ────────────────────────────────────────────
    st.write("**[4/5]** Predicting protein structure via ESMFold…")
    t4 = time.time()
    try:
        ref_pdb, mut_pdb, splice_info = _run_protein(
            ref_seq, mut_seq, win_start, edit_start, edit_end, mod_type, chrom,
            ref_out=ref_out, mut_out=mut_out,
        )
        pdb_status = (
            f"ref PDB = {len(ref_pdb):,} chars  · mut PDB = {len(mut_pdb):,} chars"
            if ref_pdb else "no coding gene in window"
        )
        splice_flag = ""
        if splice_info.get("splice_disrupted"):
            n = len(splice_info.get("skipped_exons", []))
            splice_flag = f"  ⚠ splice disruption ({n} exon(s) skipped)"
        st.write(f"  ✅ {pdb_status}{splice_flag}  · _{time.time() - t4:.1f}s_")
    except Exception as exc:
        ref_pdb = mut_pdb = ""
        splice_info = {"error": str(exc)}
        st.write(f"  ⚠ Protein folding failed: `{exc}`  · _{time.time() - t4:.1f}s_")

    # ── Step 5: Evo 2 ─────────────────────────────────────────────────────
    st.write("**[5/5]** Running Evo 2 variant scoring (k-mer log-odds)…")
    t5 = time.time()
    evo2_result, evo2_scan = _run_evo2(ref_seq, mut_seq)
    score = evo2_result.get("score", float("nan"))
    backend = evo2_result.get("backend", "kmer_fallback")
    st.write(
        f"  ✅ Score = {score:.4f}  · backend = `{backend}`  "
        f"· _{time.time() - t5:.1f}s_"
    )

    elapsed = time.time() - t0
    st.write(f"**Pipeline complete in {elapsed:.1f}s** 🎉")

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
            "evo2":      evo2_result,
            "evo2_scan": evo2_scan,
            "splice":    splice_info,
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


def _project_coords_2d(coords: np.ndarray) -> np.ndarray:
    if coords.size == 0:
        return np.zeros((0, 2))
    centered = coords - coords.mean(axis=0, keepdims=True)
    _, _, vt = np.linalg.svd(centered, full_matrices=False)
    return centered @ vt[:2].T


def _resample_track(values: np.ndarray, n_out: int) -> np.ndarray:
    if n_out <= 0:
        return np.array([])
    arr = np.asarray(values, dtype=float)
    if arr.size == 0:
        return np.zeros(n_out)
    if arr.size == n_out:
        return arr
    x_old = np.linspace(0.0, 1.0, arr.size)
    x_new = np.linspace(0.0, 1.0, n_out)
    return np.interp(x_new, x_old, arr)


def _distance_edges(coords: np.ndarray, quantile: float = 0.18,
                    min_separation: int = 2) -> list[tuple[int, int, float]]:
    if len(coords) < 3:
        return []
    pairs = []
    dists = []
    for i in range(len(coords)):
        for j in range(i + min_separation, len(coords)):
            dist = float(np.linalg.norm(coords[i] - coords[j]))
            pairs.append((i, j, dist))
            dists.append(dist)
    if not dists:
        return []
    thresh = float(np.quantile(np.array(dists), quantile))
    kept = [(i, j, d) for i, j, d in pairs if d <= thresh]
    return sorted(kept, key=lambda x: x[2])


def _map_position_to_index(position: float, win_start: int, win_end: int,
                           n_bins: int) -> int:
    if n_bins <= 1 or win_end <= win_start:
        return 0
    frac = (position - win_start) / max(1, (win_end - win_start))
    idx = int(round(frac * (n_bins - 1)))
    return max(0, min(n_bins - 1, idx))


def _make_edge_traces(node_xy: np.ndarray, edges: list[tuple[int, int, float]],
                      color: str = "rgba(180,180,180,0.35)", width_scale: float = 2.5):
    traces = []
    if not edges:
        return traces
    inv_dists = np.array([1.0 / max(d, 1e-6) for _, _, d in edges], dtype=float)
    lo, hi = float(inv_dists.min()), float(inv_dists.max())
    for (i, j, d), inv in zip(edges, inv_dists):
        norm = 0.5 if hi == lo else (inv - lo) / (hi - lo)
        traces.append(
            go.Scatter(
                x=[node_xy[i, 0], node_xy[j, 0]],
                y=[node_xy[i, 1], node_xy[j, 1]],
                mode="lines",
                line=dict(color=color, width=1.0 + norm * width_scale),
                hoverinfo="text",
                text=[f"Edge {i} ↔ {j}<br>distance={d:.3f}"] * 2,
                showlegend=False,
            )
        )
    return traces


def _render_accessibility_graph(result: dict, edit_start: int, edit_end: int):
    try:
        _, ref_coords, mut_coords = run_chromatin_sim(result["ref_hic"], result["mut_hic"])
    except Exception as exc:
        st.info(f"Accessibility graph unavailable: {exc}")
        return

    n = len(ref_coords)
    if n == 0:
        st.info("Accessibility graph unavailable for this locus.")
        return

    ref_xy = _project_coords_2d(ref_coords)
    mut_xy = _project_coords_2d(mut_coords)
    ref_atac = _resample_track(result["ref_atac"], n)
    mut_atac = _resample_track(result["mut_atac"], n)
    delta_atac = np.log2((ref_atac + 1e-6) / (mut_atac + 1e-6))
    bead_positions = np.linspace(result["win_start"], result["win_start"] + len(result["ref_seq"]), n)
    bead_in_edit = (bead_positions >= edit_start) & (bead_positions <= edit_end)
    cmin = float(delta_atac.min())
    cmax = float(delta_atac.max())
    if np.isclose(cmin, cmax):
        cmin -= 1e-6
        cmax += 1e-6

    ref_edges = _distance_edges(ref_coords)
    mut_edges = _distance_edges(mut_coords)

    def _fig_for(coords_2d: np.ndarray, edges: list[tuple[int, int, float]], title: str):
        fig = go.Figure()
        for trace in _make_edge_traces(coords_2d, edges):
            fig.add_trace(trace)
        fig.add_trace(
            go.Scatter(
                x=coords_2d[:, 0],
                y=coords_2d[:, 1],
                mode="markers+text",
                text=[str(i) if bead_in_edit[i] else "" for i in range(n)],
                textposition="top center",
                marker=dict(
                    size=[14 if bead_in_edit[i] else 10 for i in range(n)],
                    color=delta_atac,
                    colorscale="RdBu_r",
                    cmin=cmin,
                    cmax=cmax,
                    line=dict(
                        width=[2 if bead_in_edit[i] else 0.6 for i in range(n)],
                        color=["#ff5a5a" if bead_in_edit[i] else "rgba(20,20,20,0.35)" for i in range(n)],
                    ),
                    colorbar=dict(title="ATAC Δ<br>log2(Ref/Mut)"),
                ),
                hovertemplate=(
                    "Bead %{customdata[0]}<br>"
                    "Genomic pos: %{customdata[1]:,.0f}<br>"
                    "Ref ATAC: %{customdata[2]:.3f}<br>"
                    "Mut ATAC: %{customdata[3]:.3f}<br>"
                    "ATAC Δ: %{customdata[4]:+.3f}<extra></extra>"
                ),
                customdata=np.column_stack([np.arange(n), bead_positions, ref_atac, mut_atac, delta_atac]),
                showlegend=False,
            )
        )
        fig.update_layout(
            title=title,
            height=520,
            plot_bgcolor="white",
            paper_bgcolor="white",
            margin=dict(t=40, b=20, l=20, r=20),
            xaxis=dict(showgrid=False, zeroline=False, visible=False),
            yaxis=dict(showgrid=False, zeroline=False, visible=False, scaleanchor="x", scaleratio=1),
        )
        return fig

    tabs = st.tabs(["Reference accessibility graph", "Mutant accessibility graph"])
    with tabs[0]:
        st.plotly_chart(
            _fig_for(ref_xy, ref_edges, "Accessibility Graph — Reference spatial contacts"),
            use_container_width=True,
        )
    with tabs[1]:
        st.plotly_chart(
            _fig_for(mut_xy, mut_edges, "Accessibility Graph — Mutant spatial contacts"),
            use_container_width=True,
        )

    st.caption(
        "Nodes are polymer beads colored by ATAC change. Edges connect only strong spatial contacts in the simulated 3D conformation."
    )


def _build_regulatory_nodes(genes_df: pd.DataFrame, result: dict,
                            edit_start: int, edit_end: int) -> pd.DataFrame:
    if genes_df is None or genes_df.empty:
        return pd.DataFrame()

    win_start = result["win_start"]
    win_end = win_start + len(result["ref_seq"])
    n_contact = result["ref_hic"].shape[0]
    ref_atac = _resample_track(result["ref_atac"], n_contact)
    mut_atac = _resample_track(result["mut_atac"], n_contact)
    ref_cage = _resample_track(result["multi_tracks"]["CAGE"][0], n_contact)
    mut_cage = _resample_track(result["multi_tracks"]["CAGE"][1], n_contact)
    ref_ctcf = _resample_track(result["multi_tracks"]["CTCF"][0], n_contact)
    mut_ctcf = _resample_track(result["multi_tracks"]["CTCF"][1], n_contact)

    rows = []
    for _, row in genes_df.iterrows():
        raw_type = str(row["Type"])
        if raw_type == "Gene":
            node_type = "Promoter"
            pos = float(row["Start"])
        elif raw_type in {"Promoter"}:
            node_type = "Promoter"
            pos = float((row["Start"] + row["End"]) / 2)
        elif raw_type in {"Enhancer", "cCRE"}:
            node_type = "Enhancer"
            pos = float((row["Start"] + row["End"]) / 2)
        elif raw_type in {"CTCF", "Insulator"}:
            node_type = "CTCF anchor"
            pos = float((row["Start"] + row["End"]) / 2)
        else:
            continue

        idx = _map_position_to_index(pos, win_start, win_end, n_contact)
        if node_type == "Promoter":
            ref_activity = 0.7 * ref_cage[idx] + 0.3 * ref_atac[idx]
            mut_activity = 0.7 * mut_cage[idx] + 0.3 * mut_atac[idx]
        elif node_type == "Enhancer":
            ref_activity = ref_atac[idx]
            mut_activity = mut_atac[idx]
        else:
            ref_activity = 0.7 * ref_ctcf[idx] + 0.3 * ref_atac[idx]
            mut_activity = 0.7 * mut_ctcf[idx] + 0.3 * mut_atac[idx]

        rows.append(
            {
                "Name": row["Name"],
                "NodeType": node_type,
                "Position": pos,
                "Index": idx,
                "RefActivity": float(ref_activity),
                "MutActivity": float(mut_activity),
                "DeltaActivity": float(np.log2((ref_activity + 1e-6) / (mut_activity + 1e-6))),
                "OverlapsEdit": bool(pos >= edit_start and pos <= edit_end),
            }
        )

    if not rows:
        return pd.DataFrame()

    nodes = pd.DataFrame(rows)
    nodes = (
        nodes.sort_values(
            by=["NodeType", "RefActivity"],
            ascending=[True, False],
        )
        .groupby("NodeType", group_keys=False)
        .head(7)
        .reset_index(drop=True)
    )
    return nodes


def _render_regulatory_graph(genes_df: pd.DataFrame, result: dict,
                             edit_start: int, edit_end: int):
    nodes = _build_regulatory_nodes(genes_df, result, edit_start, edit_end)
    if nodes.empty:
        st.info("Regulatory graph unavailable: no suitable annotated promoters, enhancers, or CTCF anchors in view.")
        return

    _, ref_coords, mut_coords = run_chromatin_sim(result["ref_hic"], result["mut_hic"])
    n = len(ref_coords)
    ref_xy = _project_coords_2d(ref_coords)
    mut_xy = _project_coords_2d(mut_coords)
    win_start = result["win_start"]
    win_end = win_start + len(result["ref_seq"])

    node_xy_ref = np.array([ref_xy[_map_position_to_index(p, win_start, win_end, n)] for p in nodes["Position"]])
    node_xy_mut = np.array([mut_xy[_map_position_to_index(p, win_start, win_end, n)] for p in nodes["Position"]])
    node_xy = 0.5 * (node_xy_ref + node_xy_mut)

    ref_contact = result["ref_hic"]
    mut_contact = result["mut_hic"]
    edges = []
    for i in range(len(nodes)):
        for j in range(i + 1, len(nodes)):
            idx_i = int(nodes.iloc[i]["Index"])
            idx_j = int(nodes.iloc[j]["Index"])
            if idx_i == idx_j:
                continue
            ref_c = float(ref_contact[idx_i, idx_j]) if idx_i < ref_contact.shape[0] and idx_j < ref_contact.shape[1] else 0.0
            mut_c = float(mut_contact[idx_i, idx_j]) if idx_i < mut_contact.shape[0] and idx_j < mut_contact.shape[1] else 0.0
            ref_score = float(np.sqrt(max(nodes.iloc[i]["RefActivity"], 0.0) * max(nodes.iloc[j]["RefActivity"], 0.0)) * ref_c)
            mut_score = float(np.sqrt(max(nodes.iloc[i]["MutActivity"], 0.0) * max(nodes.iloc[j]["MutActivity"], 0.0)) * mut_c)
            max_score = max(ref_score, mut_score)
            if max_score <= 0:
                continue
            edges.append(
                {
                    "i": i,
                    "j": j,
                    "RefScore": ref_score,
                    "MutScore": mut_score,
                    "DeltaScore": mut_score - ref_score,
                    "MaxScore": max_score,
                }
            )

    if not edges:
        st.info("Regulatory graph unavailable: no activity-weighted contact edges passed threshold.")
        return

    edge_df = pd.DataFrame(edges)
    thresh = float(edge_df["MaxScore"].quantile(0.75))
    edge_df = edge_df[edge_df["MaxScore"] >= thresh].sort_values("MaxScore", ascending=False).head(18)
    if edge_df.empty:
        st.info("Regulatory graph unavailable: no strong rewiring edges passed threshold.")
        return

    fig = go.Figure()
    abs_delta = max(float(edge_df["DeltaScore"].abs().max()), 1e-6)
    for _, edge in edge_df.iterrows():
        i = int(edge["i"])
        j = int(edge["j"])
        delta = float(edge["DeltaScore"])
        width = 1.5 + 6.0 * (edge["MaxScore"] / max(float(edge_df["MaxScore"].max()), 1e-6))
        color = "#d1495b" if delta > 0 else "#3d7ee8"
        dash = "solid" if delta > 0 else "dot"
        fig.add_trace(
            go.Scatter(
                x=[node_xy[i, 0], node_xy[j, 0]],
                y=[node_xy[i, 1], node_xy[j, 1]],
                mode="lines",
                line=dict(color=color, width=width, dash=dash),
                hoverinfo="text",
                text=[
                    f"{nodes.iloc[i]['Name']} ↔ {nodes.iloc[j]['Name']}<br>"
                    f"Ref score: {edge['RefScore']:.3f}<br>"
                    f"Mut score: {edge['MutScore']:.3f}<br>"
                    f"Δ score: {delta:+.3f}"
                ] * 2,
                showlegend=False,
            )
        )

    type_color = {"Promoter": "#d7a630", "Enhancer": "#0f766e", "CTCF anchor": "#6b7280"}
    fig.add_trace(
        go.Scatter(
            x=node_xy[:, 0],
            y=node_xy[:, 1],
            mode="markers+text",
            text=nodes["Name"],
            textposition="top center",
            marker=dict(
                size=12 + 18 * (nodes[["RefActivity", "MutActivity"]].max(axis=1) / max(nodes[["RefActivity", "MutActivity"]].max().max(), 1e-6)),
                color=[type_color.get(t, "#0f766e") for t in nodes["NodeType"]],
                line=dict(
                    width=[2 if flag else 0.8 for flag in nodes["OverlapsEdit"]],
                    color=["#ff5a5a" if flag else "rgba(20,20,20,0.4)" for flag in nodes["OverlapsEdit"]],
                ),
            ),
            hovertemplate=(
                "%{customdata[0]}<br>"
                "Type: %{customdata[1]}<br>"
                "Ref activity: %{customdata[2]:.3f}<br>"
                "Mut activity: %{customdata[3]:.3f}<br>"
                "Δ activity: %{customdata[4]:+.3f}<extra></extra>"
            ),
            customdata=np.column_stack(
                [
                    nodes["Name"],
                    nodes["NodeType"],
                    nodes["RefActivity"],
                    nodes["MutActivity"],
                    nodes["DeltaActivity"],
                ]
            ),
            showlegend=False,
        )
    )
    fig.update_layout(
        title="Regulatory Rewiring Graph — activity × contact edges",
        height=560,
        plot_bgcolor="white",
        paper_bgcolor="white",
        margin=dict(t=50, b=20, l=20, r=20),
        xaxis=dict(showgrid=False, zeroline=False, visible=False),
        yaxis=dict(showgrid=False, zeroline=False, visible=False, scaleanchor="x", scaleratio=1),
    )
    st.plotly_chart(fig, use_container_width=True)
    st.caption(
        "Promoters, enhancers, and CTCF-like anchors are linked by activity × contact scores. Red edges are gained in the mutant; blue edges are weakened or lost."
    )


def _compute_insulation_profile(contact_map: np.ndarray, window_bins: int = 4) -> np.ndarray:
    n = contact_map.shape[0]
    profile = np.full(n, np.nan)
    for i in range(window_bins, n - window_bins):
        block = contact_map[i - window_bins:i, i:i + window_bins]
        profile[i] = float(np.nanmean(block))
    if np.all(np.isnan(profile)):
        return np.zeros(n)
    profile = pd.Series(profile).interpolate(limit_direction="both").fillna(method="bfill").fillna(method="ffill").to_numpy()
    return profile


def _estimate_boundary_consequence(result: dict, edit_start: int, edit_end: int) -> dict:
    ref_hic = np.asarray(result["ref_hic"], dtype=float)
    mut_hic = np.asarray(result["mut_hic"], dtype=float)
    n = ref_hic.shape[0]
    if n == 0:
        return {}
    edit_mid = 0.5 * (edit_start + edit_end)
    idx = _map_position_to_index(edit_mid, result["win_start"], result["win_start"] + len(result["ref_seq"]), n)
    window_bins = max(3, min(8, n // 20 if n >= 20 else 3))
    lo = max(window_bins, idx - 3)
    hi = min(n - window_bins, idx + 4)

    ref_ins = _compute_insulation_profile(ref_hic, window_bins=window_bins)
    mut_ins = _compute_insulation_profile(mut_hic, window_bins=window_bins)
    local_indices = np.arange(lo, hi)
    if len(local_indices) == 0:
        local_indices = np.array([idx])
    boundary_idx = int(local_indices[np.nanargmin(ref_ins[local_indices])])

    def _crossing_score(mat: np.ndarray, center: int) -> float:
        left0 = max(0, center - window_bins)
        left1 = center
        right0 = center
        right1 = min(mat.shape[0], center + window_bins)
        if left1 <= left0 or right1 <= right0:
            return 0.0
        return float(np.nanmean(mat[left0:left1, right0:right1]))

    ref_cross = _crossing_score(ref_hic, boundary_idx)
    mut_cross = _crossing_score(mut_hic, boundary_idx)
    ins_ref = float(ref_ins[boundary_idx])
    ins_mut = float(mut_ins[boundary_idx])
    ins_delta = ins_mut - ins_ref
    cross_delta = mut_cross - ref_cross
    if ins_delta > 0 and cross_delta > 0:
        state = "Boundary weakened"
    elif ins_delta < 0 and cross_delta < 0:
        state = "Boundary strengthened"
    else:
        state = "Boundary mixed / ambiguous"

    genomic_pos = int(np.linspace(result["win_start"], result["win_start"] + len(result["ref_seq"]), n)[boundary_idx])
    return {
        "boundary_idx": boundary_idx,
        "boundary_position": genomic_pos,
        "insulation_ref": ins_ref,
        "insulation_mut": ins_mut,
        "insulation_delta": ins_delta,
        "crossing_ref": ref_cross,
        "crossing_mut": mut_cross,
        "crossing_delta": cross_delta,
        "state": state,
    }


def _build_affected_elements_df(genes_df: pd.DataFrame, result: dict,
                                edit_start: int, edit_end: int) -> pd.DataFrame:
    if genes_df is None or genes_df.empty:
        return pd.DataFrame(columns=["Element", "Type", "Position", "ATAC Δ", "CAGE Δ", "Contact Δ", "Priority"])
    n_track = len(result["ref_atac"])
    n_contact = result["ref_hic"].shape[0]
    atac_ref = np.asarray(result["ref_atac"], dtype=float)
    atac_mut = np.asarray(result["mut_atac"], dtype=float)
    cage_ref = np.asarray(result["multi_tracks"]["CAGE"][0], dtype=float)
    cage_mut = np.asarray(result["multi_tracks"]["CAGE"][1], dtype=float)
    ref_contact = np.asarray(result["ref_hic"], dtype=float)
    mut_contact = np.asarray(result["mut_hic"], dtype=float)
    rows = []
    for _, row in genes_df.iterrows():
        mid = float((row["Start"] + row["End"]) / 2)
        i_track = _map_position_to_index(mid, result["win_start"], result["win_start"] + len(result["ref_seq"]), n_track)
        i_contact = _map_position_to_index(mid, result["win_start"], result["win_start"] + len(result["ref_seq"]), n_contact)
        atac_delta = float(np.log2((atac_ref[i_track] + 1e-6) / (atac_mut[i_track] + 1e-6)))
        cage_delta = float(np.log2((cage_ref[i_track] + 1e-6) / (cage_mut[i_track] + 1e-6)))
        contact_delta = float(np.nanmean(mut_contact[i_contact]) - np.nanmean(ref_contact[i_contact]))
        priority = abs(atac_delta) + 0.7 * abs(cage_delta) + 0.25 * abs(contact_delta)
        rows.append(
            {
                "Element": row["Name"],
                "Type": row["Type"],
                "Position": int(mid),
                "ATAC Δ": round(atac_delta, 3),
                "CAGE Δ": round(cage_delta, 3),
                "Contact Δ": round(contact_delta, 3),
                "Overlaps Edit": bool(row["Start"] <= edit_end and row["End"] >= edit_start),
                "Priority": round(priority, 3),
            }
        )
    df = pd.DataFrame(rows)
    if df.empty:
        return df
    return df.sort_values("Priority", ascending=False).reset_index(drop=True)


def _build_affected_genes_df(genes_df: pd.DataFrame, result: dict) -> pd.DataFrame:
    if genes_df is None or genes_df.empty:
        return pd.DataFrame(columns=["Gene", "Promoter ATAC Δ", "Promoter CAGE Δ", "Contact Δ", "Priority"])
    gene_df = genes_df[genes_df["Type"] == "Gene"].copy()
    if gene_df.empty:
        return pd.DataFrame(columns=["Gene", "Promoter ATAC Δ", "Promoter CAGE Δ", "Contact Δ", "Priority"])
    n_track = len(result["ref_atac"])
    n_contact = result["ref_hic"].shape[0]
    atac_ref = np.asarray(result["ref_atac"], dtype=float)
    atac_mut = np.asarray(result["mut_atac"], dtype=float)
    cage_ref = np.asarray(result["multi_tracks"]["CAGE"][0], dtype=float)
    cage_mut = np.asarray(result["multi_tracks"]["CAGE"][1], dtype=float)
    ref_contact = np.asarray(result["ref_hic"], dtype=float)
    mut_contact = np.asarray(result["mut_hic"], dtype=float)
    rows = []
    for _, row in gene_df.iterrows():
        promoter_pos = float(row["Start"])
        i_track = _map_position_to_index(promoter_pos, result["win_start"], result["win_start"] + len(result["ref_seq"]), n_track)
        i_contact = _map_position_to_index(promoter_pos, result["win_start"], result["win_start"] + len(result["ref_seq"]), n_contact)
        atac_delta = float(np.log2((atac_ref[i_track] + 1e-6) / (atac_mut[i_track] + 1e-6)))
        cage_delta = float(np.log2((cage_ref[i_track] + 1e-6) / (cage_mut[i_track] + 1e-6)))
        contact_delta = float(np.nanmean(mut_contact[i_contact]) - np.nanmean(ref_contact[i_contact]))
        priority = 1.2 * abs(cage_delta) + 0.6 * abs(atac_delta) + 0.25 * abs(contact_delta)
        rows.append(
            {
                "Gene": row["Name"],
                "Promoter ATAC Δ": round(atac_delta, 3),
                "Promoter CAGE Δ": round(cage_delta, 3),
                "Contact Δ": round(contact_delta, 3),
                "Priority": round(priority, 3),
            }
        )
    return pd.DataFrame(rows).sort_values("Priority", ascending=False).reset_index(drop=True)


def _build_step3_export_payload(result: dict, genes_df: pd.DataFrame,
                                edit_start: int, edit_end: int,
                                insight: str, boundary_info: dict,
                                elements_df: pd.DataFrame, genes_ranked_df: pd.DataFrame) -> dict:
    splice = result["extras"].get("splice", {})
    evo2_data = result["extras"].get("evo2", {})
    return {
        "input": {
            "gene": st.session_state.gene,
            "chrom": st.session_state.chrom,
            "edit_start": edit_start,
            "edit_end": edit_end,
            "modification": st.session_state.mod_type,
        },
        "summary": {
            "mechanistic_insight": insight,
            "delta_stats": result["delta_stats"],
            "boundary": boundary_info,
            "protein_target": splice.get("gene"),
            "splice": splice,
            "evo2": evo2_data,
        },
        "top_elements": elements_df.head(15).to_dict(orient="records"),
        "top_genes": genes_ranked_df.head(15).to_dict(orient="records"),
        "annotations": genes_df.to_dict(orient="records") if genes_df is not None else [],
    }


# ---------------------------------------------------------------------------
# Helper: render differentiation trajectory
# ---------------------------------------------------------------------------
def _render_trajectory(result: dict):
    with st.status("Running differentiation trajectory (4 AlphaGenome calls)…",
                   expanded=True) as _traj_status:
        try:
            maps = _run_trajectory(result["ref_seq"], result["mut_seq"])
            _traj_status.update(label="Trajectory complete!", state="complete",
                                expanded=False)
        except Exception as exc:
            _traj_status.update(label=f"Trajectory failed: {exc}", state="error",
                                expanded=True)
            st.error(f"Differentiation trajectory failed: {exc}")
            return

    titles = ["H1-hESC  Ref", "H1-hESC  Mut", "GM12878  Ref", "GM12878  Mut"]
    keys   = ["h1_ref", "h1_mut", "gm12878_ref", "gm12878_mut"]
    vmax   = max(np.nanmax(np.abs(maps[k])) for k in keys if maps[k].size > 0)

    fig, axes = plt.subplots(1, 4, figsize=(16, 4))
    for ax, key, title in zip(axes, keys, titles):
        im = ax.imshow(maps[key], cmap="YlOrRd", vmin=0, vmax=vmax, aspect="auto")
        ax.set_title(title, fontsize=11)
        ax.axis("off")
    fig.colorbar(im, ax=axes.tolist(), shrink=0.6, label="Contact frequency")
    fig.text(0.5, -0.02,
             "← Embryonic stem cell (H1-hESC)   ·   Lymphoblastoid (GM12878) →",
             ha="center", fontsize=10, color="#444")
    plt.tight_layout()
    with tempfile.NamedTemporaryFile(suffix=".png", delete=False) as tmp:
        fig.savefig(tmp.name, dpi=150, bbox_inches="tight")
        plt.close(fig)
        st.image(Image.open(tmp.name), use_container_width=True)

    with st.expander("Model & Methodology"):
        st.markdown("""
**Model:** AlphaGenome (DeepMind) with cell-type-specific ontology terms
**Cell types:** H1 (`EFO:0003042`) = H1 embryonic stem cell line; GM12878 (`EFO:0002784`) = committed lymphoblastoid (B-cell lineage)
**Note:** AlphaGenome contact maps are only available for: H1, GM12878, HFFc6, HCT116, IMR-90, HeLa, HepG2, H9, KBM-7.
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

    with st.expander("What happens in Step 1", expanded=False):
        st.markdown("""
This step defines the **genomic edit scenario** that the rest of the app will analyze.

What you choose here determines:
- which chromosome and coordinate range the app studies
- what kind of DNA change is applied
- which reference and mutant sequences are built later for prediction

The app can start from either:
- a preset **gene symbol** from the built-in registry, or
- a manually entered genomic interval like `chr1:47210000-47260000`

Then you specify the edit itself:
- **Deletion**: remove the bases between start and end
- **Insertion**: add a new sequence at the chosen site
- **SNP**: change one base to another base
- **None**: analyze the region without changing sequence

Nothing biological is predicted yet in this step. This page is just building the exact DNA-edit request that will drive the later AlphaGenome, RNA, protein, and 3D analyses.
""")

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
    insert_buf = st.session_state.insert_seq
    snp_buf = st.session_state.snp_alt
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
        if mod_type == "Insertion":
            insert_buf = st.text_area(
                f"Inserted sequence ({edit_size} bp)",
                value=insert_buf,
                height=70,
                help=f"Exactly {edit_size} bases (A/C/G/T), or leave empty to use A×{edit_size}.",
            )
        if mod_type == "SNP":
            _bases = ["A", "C", "G", "T"]
            _si = _bases.index(snp_buf) if snp_buf in _bases else 0
            snp_buf = st.selectbox("Alternate base (SNP)", _bases, index=_si)
    elif edit_size <= 0:
        st.error("Edit end must be greater than edit start.")

    with st.expander("How the edit definition is used later", expanded=False):
        st.markdown("""
Once you click **Preview Locus**, the app stores this edit definition in Streamlit session state.

Later, in the prediction step, the code will:
1. fetch the real reference DNA sequence from the selected genomic window
2. apply the edit you defined here to create a mutant DNA sequence
3. run both reference and mutant sequences through the analysis pipeline

So this page is effectively the **input form for building the two DNA worlds** the app compares:
- **reference** genome sequence
- **mutant** edited sequence
""")

    st.write("")
    if st.button("Preview Locus →", type="primary", disabled=(edit_size <= 0)):
        st.session_state.chrom = chrom
        st.session_state.edit_start = int(edit_start)
        st.session_state.edit_end = int(edit_end)
        st.session_state.mod_type = mod_type
        st.session_state.insert_seq = insert_buf if mod_type == "Insertion" else ""
        st.session_state.snp_alt = snp_buf if mod_type == "SNP" else st.session_state.snp_alt
        st.session_state.locus_annotations = None   # force refresh
        st.session_state.step = 2
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

    with st.expander("What happens in Step 2", expanded=False):
        st.markdown("""
This step is the **locus preview and annotation sanity check**.

Before running heavy prediction models, the app asks:
- what annotated genes and regulatory elements are already known here?
- does the planned edit overlap anything important?
- is the edit sitting in a promoter, enhancer-rich region, gene body, or a relatively empty interval?

This helps you inspect the local genomic context before running the full HG-DT pipeline.
""")

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

    with st.expander("How the preview browser works", expanded=False):
        st.markdown("""
The genome preview is built from **real locus annotations**, not model predictions.

At this stage, the app fetches annotations from UCSC and shows:
- genes
- other annotated elements returned by the locus-annotation helper
- the planned edit region highlighted in context

The purpose is to show where your edit sits relative to known biology before prediction begins.
""")

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

    with st.expander("Why overlap checking matters", expanded=False):
        st.markdown("""
This overlap table is a quick interpretability check.

If the edit overlaps an annotated element, that gives an immediate first guess about mechanism:
- overlap with a **gene** may suggest coding or transcript effects
- overlap with a **promoter / enhancer / cCRE / CTCF-like element** may suggest regulatory or structural effects
- no obvious overlap does **not** mean the edit is irrelevant; it may still change sequence features or long-range contacts

The app uses this page to help you decide whether the edit looks biologically sensible before spending time on the full model run.
""")

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
        with st.status("Running HG-DT pipeline…", expanded=True) as _pipeline_status:
            try:
                st.session_state.pipeline_result = _run_pipeline()
                _pipeline_status.update(
                    label="Pipeline complete!",
                    state="complete",
                    expanded=False,
                )
            except Exception as exc:
                _pipeline_status.update(
                    label=f"Pipeline failed: {exc}",
                    state="error",
                    expanded=True,
                )
                st.error(f"Pipeline error: {exc}")
                st.stop()

    result   = st.session_state.pipeline_result
    _annot   = st.session_state.locus_annotations
    genes_df = _annot if _annot is not None else pd.DataFrame(
        columns=["Name", "Type", "Start", "End"]
    )
    boundary_info = _estimate_boundary_consequence(result, edit_start, edit_end)
    elements_ranked_df = _build_affected_elements_df(genes_df, result, edit_start, edit_end)
    genes_ranked_df = _build_affected_genes_df(genes_df, result)

    # ── Summary header ─────────────────────────────────────────────────────
    st.subheader("Most Likely Mechanism")
    ds = result["delta_stats"]
    mod_details = {
        "type":   st.session_state.mod_type.lower(),
        "target": f"regulatory element near {st.session_state.gene}",
    }
    insight = generate_mechanistic_insight(mod_details, ds)
    st.success(insight)

    top_gene = genes_ranked_df.iloc[0]["Gene"] if not genes_ranked_df.empty else "N/A"
    top_element = elements_ranked_df.iloc[0]["Element"] if not elements_ranked_df.empty else "N/A"
    coding_note = result["extras"].get("splice", {}).get("gene") or "No coding target"

    s1, s2, s3, s4 = st.columns(4)
    with s1:
        st.metric("Top affected gene", str(top_gene))
    with s2:
        st.metric("Top affected element", str(top_element))
    with s3:
        st.metric("Boundary state", boundary_info.get("state", "N/A"))
    with s4:
        st.metric("Protein target", str(coding_note))

    with st.expander("How this summary is assembled", expanded=False):
        st.markdown("""
This header combines the strongest signals already computed in the app:

- accessibility and expression deltas
- loop weakening / strengthening
- boundary consequence near the edit
- top-ranked affected genes and annotated elements
- the coding target chosen for the protein branch

It is meant to give you a fast answer to: **what is the most likely biological story here?**
""")

    # ── Two-column layout ──────────────────────────────────────────────────
    col_sys1, col_sys2 = st.columns(2, gap="large")

    with col_sys1:
        st.subheader("System 1 — Protein Path")
        splice = result["extras"].get("splice", {})
        gene_name = splice.get("gene", "")
        ref_len   = splice.get("ref_aa_len", 0)
        mut_len   = splice.get("mut_aa_len", 0)
        if gene_name:
            st.caption(
                f"Gene: **{gene_name}** · "
                f"Ref: {ref_len} AA · Mut: {mut_len} AA · "
                f"Δ {mut_len - ref_len:+d} AA"
            )
        else:
            st.caption("DNA → mRNA isoforms → Amino acid translation → 3D structure")

        fc = splice.get("fold_cap_aa") or 400
        if splice.get("truncated_to_fold"):
            st.warning(
                f"ESMFold input is capped at **{fc}** residues (Meta public API limit). "
                f"Longer chains are truncated; full-length folding (e.g. AlphaFold 3) is planned."
            )

        with st.expander("How the protein prediction works", expanded=False):
            st.markdown(f"""
This panel follows a **DNA -> RNA -> protein -> 3D structure** path for the coding gene nearest the edit.

1. **Find nearby gene models:** the app fetches real UCSC RefSeq gene models in the prediction window.
2. **Choose a coding target:** it picks a coding transcript that overlaps the edit if possible, otherwise the nearest coding gene.
3. **Assemble transcript sequence:** exon segments are stitched together from the reference and mutant DNA windows.
4. **Check splice disruption:** AlphaGenome splice-site scores are used to flag exons whose donor/acceptor sites may be weakened enough to skip.
5. **Convert transcript DNA to mRNA:** the exon-only transcript is strand-corrected and `T` bases are converted to `U`.
6. **Find the longest ORF:** the app scans the mRNA for the longest start-to-stop coding region.
7. **Translate to amino acids:** that ORF is translated with Biopython into a protein sequence.
8. **Fold the protein:** the amino acid sequence is sent to ESMFold, then the reference and mutant structures are overlaid here.

In this run, the protein-focused target is **{gene_name or 'the nearest coding gene'}**.
""")

        render_protein_overlay(result["ref_pdb"], result["mut_pdb"])

        if splice.get("splice_disrupted"):
            n_skipped = len(splice.get("skipped_exons", []))
            st.warning(
                f"**{gene_name}**: length change detected — "
                f"ref {ref_len} AA → mut {mut_len} AA"
                + (", frameshift" if splice.get("frameshift") else "") + "."
            )
        elif splice and not splice.get("error"):
            overlap_note = ""
            edit_mid = (edit_start + edit_end) / 2
            # Warn if the chosen gene doesn't actually overlap the edit
            if gene_name and gene_name != st.session_state.gene:
                overlap_note = (
                    f" (nearest coding gene — **{st.session_state.gene}** "
                    f"may be non-coding or fully deleted)"
                )
            st.success(f"No coding sequence disruption in **{gene_name}**." + overlap_note)
        elif splice.get("error"):
            st.info(f"Protein prediction note: {splice['error']}")

        with st.expander("Model & Methodology"):
            st.markdown("""
**Structure:** NVIDIA ESMFold NIM → Meta public ESMFold (free), max ~400 AA
**RNA:** Longest-ORF translation from exon-assembled mRNA (`predict_isoforms` + `compare_translation`), with CDS fallback
**Overlay:** Blue = Reference · Orange = Mutant (both in one py3Dmol scene)
**Splicing:** AlphaGenome SPLICE_SITES track → donor/acceptor score per exon
**Gene models:** UCSC hg38 RefSeq (live REST API)
**Confidence:** B-factor column = pLDDT (0–100). Regions <50 = unreliably folded.
""")

    with col_sys2:
        st.subheader("System 2 — Regulatory Path")
        st.caption("DNA → Chromatin accessibility → Enhancer loops → Expression")

        with st.expander("How the regulatory prediction works", expanded=False):
            st.markdown("""
This panel shows the **regulatory consequence** branch of the pipeline.

1. The app builds matched **reference** and **mutant** DNA windows around the edit.
2. AlphaGenome predicts regulatory tracks on both windows:
   - **ATAC** for accessibility
   - **CAGE** as a promoter / transcription proxy
   - **CHIP_TF / CTCF** for structural binding signal
   - **CONTACT_MAPS** for 3D chromatin contacts
   - **SPLICE_SITES** for splice disruption
3. The app computes deltas between reference and mutant outputs.
4. Those deltas are summarized into accessibility loss, expression change, and loop rewiring metrics.

The track panels below are therefore not measured wet-lab signals from this exact edit; they are **model predictions** comparing the edited sequence to the original sequence.
""")

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
    with st.expander("What this accessibility plot means", expanded=False):
        st.markdown("""
This chart maps predicted accessibility change back onto **real annotated genomic elements** in the locus.

- Each bubble is an annotated element such as a gene, enhancer, promoter, CTCF site, or cCRE.
- The app samples the predicted ATAC signal around that element in the reference and mutant runs.
- The Y-axis is `log2(Ref / Mut)`, so positive values mean the element is predicted to become **less accessible** after the edit.
- Larger bubbles represent wider elements.

This is useful for seeing **which annotated elements are most affected**, rather than only looking at the raw track lines.
""")
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
        with st.expander("How to read the Hi-C contact map", expanded=False):
            st.markdown("""
This is a **difference contact map** built from AlphaGenome-predicted contact maps.

- The reference and mutant 3D contact maps are predicted separately.
- The app subtracts them to show **where contacts increase or decrease** after the edit.
- **Red** regions indicate gained contacts in the mutant.
- **Blue** regions indicate lost contacts in the mutant.
- Because the map is shown in triangular form, boundary shifts and domain rewiring are easier to spot visually.

This panel answers the question: **did the edit change 3D genome organization?**
""")

    st.subheader("3D Chromatin Conformation — Reference → Mutant")
    st.caption(
        "Polymer simulation (Langevin dynamics) on AlphaGenome contact maps. "
        "Press ▶ to watch chromatin reorganize. "
        "Rainbow = genomic position (5′→3′). Red beads = edit region."
    )
    with st.expander("How the 3D chromatin animation is generated", expanded=False):
        st.markdown("""
This animation is a **coarse-grained polymer simulation**, not an atom-by-atom chromosome model.

1. The app takes the predicted reference and mutant contact maps.
2. Each map is downsampled to a smaller number of bins so it can be simulated efficiently.
3. Strong contacts are converted into **harmonic spring restraints** between beads.
4. A bead-spring polymer is simulated with:
   - backbone springs between neighboring beads
   - contact-derived springs between distant beads
   - soft repulsion so beads do not collapse onto each other
   - thermal noise through overdamped Langevin dynamics
5. The reference and mutant 3D conformations are aligned, and the app interpolates between them into an animation.

So the movie is best interpreted as **how the locus may reorganize in 3D**, not as an exact physical reconstruction of chromatin.
""")
    window_pad = 200_000
    render_chromatin_animation(
        result["ref_hic"], result["mut_hic"],
        edit_start=edit_start,
        edit_end=edit_end,
        genomic_start=max(0, edit_start - window_pad),
        genomic_end=edit_end + window_pad,
    )

    # ── Graph views from polymer + activity ───────────────────────────────
    st.divider()
    st.subheader("Graph Views of the Locus")
    st.caption(
        "These network views reuse the polymer simulation and AlphaGenome tracks to show which regions are brought together in 3D and which annotated elements are likely rewired."
    )
    gtab1, gtab2 = st.tabs(["Accessibility graph", "Regulatory graph"])
    with gtab1:
        with st.expander("What this graph shows", expanded=False):
            st.markdown("""
This is an **accessibility graph** built on top of the polymer simulation.

- **Nodes:** polymer beads representing genomic bins
- **Node color:** predicted ATAC change, shown as `log2(Ref / Mut)`
- **Edges:** only strong spatial contacts in the simulated 3D structure

The idea is to show **which accessible or accessibility-changing regions are physically brought together** in the reference and mutant conformations.
""")
        _render_accessibility_graph(result, edit_start, edit_end)

    with gtab2:
        with st.expander("What this graph shows", expanded=False):
            st.markdown("""
This is a **regulatory rewiring graph** built from annotated genomic elements.

- **Nodes:** promoters, enhancers, and CTCF-like anchors from the locus annotations
- **Node size:** element activity strength
- **Edges:** activity × contact score
- **Edge color:** red = gained in mutant, blue = weakened or lost

This view is meant to highlight **which functional elements are most likely to be rewired** after the edit, rather than only showing raw geometry.
""")
        _render_regulatory_graph(genes_df, result, edit_start, edit_end)

    # ── Ranked tables ──────────────────────────────────────────────────────
    st.divider()
    st.subheader("Top Affected Genes and Elements")
    with st.expander("Why these tables are useful", expanded=False):
        st.markdown("""
These tables turn the model outputs into ranked candidates you can inspect or export.

- **Top affected genes** rank genes by promoter accessibility, promoter expression proxy, and local contact change.
- **Top affected elements** rank annotated genes and regulatory elements by accessibility, expression-like change, and contact rewiring.

These are useful when you want to move from visual interpretation to **experimental prioritization**.
""")
    t1, t2 = st.tabs(["Top affected genes", "Top affected elements"])
    with t1:
        if not genes_ranked_df.empty:
            st.dataframe(genes_ranked_df.head(15), use_container_width=True, hide_index=True)
        else:
            st.info("No gene-ranked table available for this locus.")
    with t2:
        if not elements_ranked_df.empty:
            st.dataframe(elements_ranked_df.head(20), use_container_width=True, hide_index=True)
        else:
            st.info("No element-ranked table available for this locus.")

    # ── Boundary consequence ───────────────────────────────────────────────
    st.divider()
    st.subheader("TAD Boundary Consequence")
    with st.expander("How this boundary diagnostic works", expanded=False):
        st.markdown("""
This card looks near the edit for the strongest local insulation minimum in the predicted contact map.

It then compares:
- **local insulation** in the reference and mutant maps
- **cross-boundary contact** across the same neighborhood

If insulation rises and cross-boundary contact also rises, that is consistent with a **weakened boundary**.
""")
    if boundary_info:
        b1, b2, b3, b4 = st.columns(4)
        with b1:
            st.metric("Boundary state", boundary_info["state"])
        with b2:
            st.metric("Boundary position", f"{boundary_info['boundary_position']:,}")
        with b3:
            st.metric("Insulation Δ", f"{boundary_info['insulation_delta']:+.3f}",
                      delta="weaker" if boundary_info["insulation_delta"] > 0 else "stronger",
                      delta_color="inverse")
        with b4:
            st.metric("Cross-boundary contact Δ", f"{boundary_info['crossing_delta']:+.3f}",
                      delta="more crossing" if boundary_info["crossing_delta"] > 0 else "less crossing",
                      delta_color="inverse")
    else:
        st.info("Boundary consequence could not be estimated for this run.")

    # ── Differentiation Trajectory ─────────────────────────────────────────
    st.divider()
    st.subheader("Differentiation Trajectory")
    st.caption(
        "AlphaGenome contact maps at H1-hESC (embryonic stem) → GM12878 (committed lymphoblastoid). "
        "Shows how the edit reshapes 3D organization across the developmental axis."
    )
    with st.expander("What the differentiation trajectory is doing", expanded=False):
        st.markdown("""
This section reruns the contact-map prediction in different supported AlphaGenome cell contexts.

The goal is to ask: **does the same sequence edit have different 3D consequences in different cell states?**

In other words, this is a quick way to compare whether the edit looks more disruptive in a stem-like context versus a more differentiated context.
""")
    _render_trajectory(result)

    # ── Causal Summary + exports ───────────────────────────────────────────
    st.divider()
    st.subheader("Causal Summary")
    with st.expander("How the summary is generated", expanded=False):
        st.markdown("""
The causal summary is a text interpretation built from a few core delta statistics:

- predicted accessibility change
- predicted expression change
- predicted loop weakening or strengthening

It is meant to summarize the most likely mechanistic story of the edit in plain language, based on the model outputs above.
""")
    st.success(insight)

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

    export_payload = _build_step3_export_payload(
        result, genes_df, edit_start, edit_end,
        insight, boundary_info, elements_ranked_df, genes_ranked_df
    )
    export_json = json.dumps(export_payload, indent=2)
    ex1, ex2, ex3 = st.columns(3)
    with ex1:
        st.download_button(
            "Download summary JSON",
            export_json,
            file_name=f"{st.session_state.gene}_{st.session_state.mod_type.lower()}_summary.json",
            mime="application/json",
            use_container_width=True,
        )
    with ex2:
        st.download_button(
            "Download genes CSV",
            genes_ranked_df.to_csv(index=False),
            file_name=f"{st.session_state.gene}_{st.session_state.mod_type.lower()}_genes.csv",
            mime="text/csv",
            use_container_width=True,
        )
    with ex3:
        st.download_button(
            "Download elements CSV",
            elements_ranked_df.to_csv(index=False),
            file_name=f"{st.session_state.gene}_{st.session_state.mod_type.lower()}_elements.csv",
            mime="text/csv",
            use_container_width=True,
        )

    ct_note = result["extras"].get("cell_type")
    if ct_note:
        st.caption(f"Predictions filtered to cell type: **{ct_note}**")
    st.caption(
        "Mechanistic attribution from AlphaGenome track deltas (log₂ Ref/Mut). "
        "Sequence from UCSC hg38. No hardcoded or mock data."
    )

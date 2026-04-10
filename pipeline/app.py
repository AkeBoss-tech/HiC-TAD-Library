#!/usr/bin/env python3
"""
Streamlit UI for the Drug Discovery Pipeline.

Run with:
    streamlit run pipeline/app.py

Separate from synthesis/app.py — has its own st.set_page_config().
"""

import importlib
import json
import os
import sys

import numpy as np
import streamlit as st

# ── Page config ────────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Drug Discovery Pipeline — UNC5B",
    page_icon="💊",
    layout="wide",
)

# Add repo root to path so pipeline.* imports work when run from any directory
_repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _repo_root not in sys.path:
    sys.path.insert(0, _repo_root)

# ── Constants ──────────────────────────────────────────────────────────────────
OUTPUTS = {
    "fasta":      "data/processed/pipeline_target_sequence.fasta",
    "embeddings": "data/processed/pipeline_esm2_embeddings.npy",
    "pdb":        "data/processed/pipeline_structure.pdb",
    "molecules":  "data/processed/pipeline_molecules.json",
    "img_s1":     "media/pipeline_stage1_expression.png",
    "img_s2":     "media/pipeline_stage2_embeddings.png",
    "img_s3":     "media/pipeline_stage3_structure.png",
    "img_s4":     "media/pipeline_stage4_molecules.png",
    "img_s4_2d":  "media/pipeline_stage4_molecules_2d.png",
}

STAGE_MODULES = {
    1: "pipeline.stage1_genomics",
    2: "pipeline.stage2_protein",
    3: "pipeline.stage3_folding",
    4: "pipeline.stage4_molecules",
}


# ── Sidebar ────────────────────────────────────────────────────────────────────

st.sidebar.title("💊 Drug Discovery Pipeline")
st.sidebar.markdown("**Target:** UNC5B death domain (Q8K1S3)")
st.sidebar.markdown("**Locus:** mm10 chr10:60.24–61.29 Mb")
st.sidebar.markdown("---")
page = st.sidebar.radio(
    "Navigate to",
    ["Overview", "Stage 1 — Genomics", "Stage 2 — Embeddings",
     "Stage 3 — Structure", "Stage 4 — Molecules", "Run Pipeline"],
)
st.sidebar.markdown("---")
st.sidebar.markdown("**Outputs**")
nvidia_key = os.getenv("NVIDIA_API_KEY") or ""
st.sidebar.metric("NVIDIA_API_KEY", "✅ Set" if nvidia_key else "❌ Missing")
for key, path in OUTPUTS.items():
    label = os.path.basename(path)
    st.sidebar.metric(label[:28], "✅" if os.path.exists(path) else "❌")


# ── Pages ──────────────────────────────────────────────────────────────────────

if page == "Overview":
    st.title("💊 In-silico Drug Discovery Pipeline")
    st.markdown("""
Automated 4-stage pipeline: **chromatin locus → protein → structure → molecules**

| Stage | Name | Tool | API calls |
|-------|------|------|-----------|
| 1 | Genomics context + sequence fetch | UniProt REST | 0 (public) |
| 2 | ESM-2 protein embeddings | NVIDIA NIM ESM2-30M | 1 |
| 3 | ESMFold structure prediction | NVIDIA NIM ESMFold | 1 |
| 4 | MolMIM molecule generation + DiffDock docking | NVIDIA NIM | 22 |

**Total:** 24 API calls · ~12 min end-to-end
""")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Target", "UNC5B death domain")
    col2.metric("Accession", "Q8K1S3 (mouse)")
    col3.metric("Domain size", "~79 aa (res 865–943)")
    col4.metric("Total API calls", "24")

    if os.path.exists(OUTPUTS["img_s1"]):
        st.image(OUTPUTS["img_s1"], caption="Stage 1 — Unc5b locus CTCF context", use_container_width=True)


elif page == "Stage 1 — Genomics":
    st.title("Stage 1 — Locus Context & Target Sequence")
    st.markdown("""
Fetches the **UNC5B death domain** sequence from UniProt (accession Q8BYU4, mouse) and
plots the AlphaGenome CTCF signal over the Unc5b locus with the gene body highlighted.

- **No API key required** — uses the public UniProt REST API and existing WT cache.
- Output: `data/processed/pipeline_target_sequence.fasta`
""")

    if os.path.exists(OUTPUTS["img_s1"]):
        st.image(OUTPUTS["img_s1"], use_container_width=True)
    else:
        st.warning("Figure not yet generated. Run **Stage 1** from the Run Pipeline page.")

    if os.path.exists(OUTPUTS["fasta"]):
        with open(OUTPUTS["fasta"]) as f:
            fasta_text = f.read()
        st.subheader("Target sequence (FASTA)")
        st.code(fasta_text, language="text")


elif page == "Stage 2 — Embeddings":
    st.title("Stage 2 — ESM-2 Protein Embeddings")
    st.markdown("""
Encodes the death domain with **ESM-2 30M** (NVIDIA NIM).
Mean-pooled embeddings (480-dim) capture evolutionary context used by MolMIM in Stage 4.
""")

    if os.path.exists(OUTPUTS["img_s2"]):
        st.image(OUTPUTS["img_s2"], use_container_width=True)
    else:
        st.warning("Figure not yet generated. Run **Stage 2** from the Run Pipeline page.")

    if os.path.exists(OUTPUTS["embeddings"]):
        emb = np.load(OUTPUTS["embeddings"])
        col1, col2 = st.columns(2)
        col1.metric("Shape", f"{emb.shape[0]} × {emb.shape[1]}")
        col2.metric("Mean |embedding|", f"{np.abs(emb).mean():.4f}")


elif page == "Stage 3 — Structure":
    st.title("Stage 3 — ESMFold Structure Prediction")
    st.markdown("""
Predicts the 3D structure of the UNC5B death domain using **ESMFold** (NVIDIA NIM).
The predicted PDB is used as the docking receptor in Stage 4.
""")

    if os.path.exists(OUTPUTS["img_s3"]):
        st.image(OUTPUTS["img_s3"], use_container_width=True)
    else:
        st.warning("Figure not yet generated. Run **Stage 3** from the Run Pipeline page.")

    if os.path.exists(OUTPUTS["pdb"]):
        with open(OUTPUTS["pdb"]) as f:
            pdb_text = f.read()
        n_atoms = sum(1 for ln in pdb_text.splitlines() if ln.startswith("ATOM"))
        n_ca    = sum(1 for ln in pdb_text.splitlines()
                      if ln.startswith("ATOM") and ln[12:16].strip() == "CA")
        col1, col2 = st.columns(2)
        col1.metric("ATOM records", n_atoms)
        col2.metric("Cα residues", n_ca)
        with st.expander("View PDB (first 30 lines)"):
            st.code("\n".join(pdb_text.splitlines()[:30]), language="text")


elif page == "Stage 4 — Molecules":
    st.title("Stage 4 — Drug Candidates")
    st.markdown("""
**MolMIM** generates 20 candidate molecules conditioned on the death domain sequence.
**DiffDock** docks each against the ESMFold structure and returns a confidence score.
Ranked by DiffDock confidence (higher = better predicted binding pose).
""")

    if os.path.exists(OUTPUTS["img_s4"]):
        st.image(OUTPUTS["img_s4"], use_container_width=True)
    else:
        st.warning("Figure not yet generated. Run **Stage 4** from the Run Pipeline page.")

    if os.path.exists(OUTPUTS["img_s4_2d"]):
        st.image(OUTPUTS["img_s4_2d"], caption="2D structures (rdkit)", use_container_width=True)

    if os.path.exists(OUTPUTS["molecules"]):
        with open(OUTPUTS["molecules"]) as f:
            molecules = json.load(f)
        st.subheader(f"Top {min(10, len(molecules))} Candidates")
        import pandas as pd
        df = pd.DataFrame(molecules[:10])
        if "smiles" in df.columns:
            df["smiles_short"] = df["smiles"].str[:40] + "..."
            df = df[["rank", "smiles_short", "molmim_score", "diffdock_confidence"]]
        st.dataframe(df, use_container_width=True)


elif page == "Run Pipeline":
    st.title("Run Pipeline")
    st.info(
        "Stages 2–4 require `NVIDIA_API_KEY` set in your `.env` file. "
        "Stage 1 uses the public UniProt REST API (no key needed)."
    )

    from dotenv import load_dotenv
    load_dotenv()

    def run_stage_ui(stage_id: int, label: str, api_calls: str, time_est: str) -> None:
        st.subheader(f"Stage {stage_id} — {label}")
        st.caption(f"API calls: {api_calls} · Estimated time: {time_est}")
        if st.button(f"▶ Run Stage {stage_id}"):
            with st.spinner(f"Running Stage {stage_id} ..."):
                try:
                    mod = importlib.import_module(STAGE_MODULES[stage_id])
                    importlib.reload(mod)
                    mod.main()
                    st.success(f"Stage {stage_id} complete!")
                    # Show generated figure if available
                    img_key = f"img_s{stage_id}"
                    if img_key in OUTPUTS and os.path.exists(OUTPUTS[img_key]):
                        st.image(OUTPUTS[img_key], use_container_width=True)
                except Exception as exc:
                    st.error(f"Stage {stage_id} failed: {exc}")

    run_stage_ui(1, "Genomics Context + Sequence Fetch",
                 "0 (UniProt)", "~5 s")
    run_stage_ui(2, "ESM-2 Protein Embeddings",
                 "1 (NVIDIA NIM ESM2-30M)", "~20 s")
    run_stage_ui(3, "ESMFold Structure Prediction",
                 "1 (NVIDIA NIM ESMFold)", "~60 s")
    run_stage_ui(4, "Molecule Generation + Docking",
                 "22 (MolMIM + 20× DiffDock)", "~10 min")

    st.subheader("Generate Report")
    if st.button("📄 Generate pipeline_report.html"):
        try:
            import pipeline.pipeline_report as pr
            importlib.reload(pr)
            pr.main()
            st.success("pipeline_report.html generated!")
            with open("pipeline_report.html") as f:
                st.download_button(
                    "⬇ Download pipeline_report.html",
                    f.read(),
                    file_name="pipeline_report.html",
                    mime="text/html",
                )
        except Exception as exc:
            st.error(f"Failed: {exc}")

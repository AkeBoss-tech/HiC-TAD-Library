#!/usr/bin/env python3
"""
Streamlit UI for Synthesis Experiments.

Run with:
    streamlit run synthesis/app.py

Lets you browse outputs, trigger individual experiments, and inspect metrics
without touching the command line for each step.
"""

import os
import sys
import warnings
import numpy as np
import streamlit as st
import matplotlib.pyplot as plt

# ── Page config ───────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="Synthesis Experiments — Chromatin Engineering",
    page_icon="⚗️",
    layout="wide",
)

# ── Constants ─────────────────────────────────────────────────────────────────
WT_CACHE  = "data/processed/synthesis_unc5b_wt_cache.npz"
MEDIA     = {
    "exp1":       "media/synthesis_exp1_baseline_unc5b.png",
    "exp2_panels":"media/synthesis_exp2_contact_panels.png",
    "exp2_delta": "media/synthesis_exp2_delta_insulation.png",
    "exp3":       "media/synthesis_exp3_synthetic_insulator.png",
}
CTCF_FWD  = "CCGCGAGGCGCAGG"
CTCF_REV  = "CCTGCGCCTCGCGG"
WIN_START = 60_240_000
WIN_END   = 61_288_576
CHROM     = "chr10"
CELL_TYPE = "EFO:0004038"


# ── Helpers ───────────────────────────────────────────────────────────────────

def load_wt_cache():
    if not os.path.exists(WT_CACHE):
        return None
    try:
        data = np.load(WT_CACHE)
        return dict(
            contact_mat=data["contact_mat"],
            ctcf_signal=data["ctcf_signal"],
            resolution=float(data["resolution"]),
            n_bins=int(data["n_bins"]),
        )
    except Exception as exc:
        st.warning(f"WT cache could not be loaded: {exc}")
        return None


def compute_insulation(mat, resolution, window_bp=50_000):
    n   = mat.shape[0]
    win = max(2, int(window_bp / resolution))
    ins = np.full(n, np.nan)
    for i in range(win, n - win):
        ins[i] = mat[i - win:i, i:i + win].mean()
    return ins


def load_model():
    from dotenv import load_dotenv
    load_dotenv()
    from alphagenome.models import dna_client
    api_key = os.getenv("ALPHA_GENOME_API_KEY")
    if not api_key:
        st.error("ALPHA_GENOME_API_KEY not set in .env — cannot call API.")
        return None, None
    return dna_client.create(api_key), dna_client


def build_ctcf_pair_seq(n_pairs, spacer=10):
    pair = CTCF_FWD + "N" * spacer + CTCF_REV
    return ("N" * spacer).join([pair] * n_pairs)


# ── Sidebar ───────────────────────────────────────────────────────────────────

st.sidebar.title("⚗️ Synthesis Experiments")
st.sidebar.markdown("**Locus:** Unc5b · chr10:60.24–61.29 Mb")
st.sidebar.markdown("**Cell type:** Mouse ESC (EFO:0004038)")
st.sidebar.markdown("---")
page = st.sidebar.radio(
    "Navigate to",
    ["Overview", "Exp 1 — Baseline", "Exp 2 — Deletions", "Exp 3 — Synthetic Insulator", "Run Experiments"],
)
st.sidebar.markdown("---")
cache_exists = os.path.exists(WT_CACHE)
st.sidebar.metric("WT cache", "✅ Present" if cache_exists else "❌ Missing")
for k, path in MEDIA.items():
    label = os.path.basename(path).replace("synthesis_", "").replace(".png", "")
    st.sidebar.metric(label, "✅" if os.path.exists(path) else "❌")


# ── Pages ─────────────────────────────────────────────────────────────────────

if page == "Overview":
    st.title("⚗️ In-silico Chromatin Boundary Engineering")
    st.markdown("""
Three validation experiments use **AlphaGenome** to probe and redesign the
*Unc5b* CTCF insulator cluster in Mouse ESC.

| Exp | What | API calls |
|-----|------|-----------|
| 1 | Wild-type baseline — contact map, CTCF signal, insulation | 1 |
| 2 | CTCF deletion series — quantify insulation loss | 3 |
| 3 | Synthetic insulator — convergent CTCF pair insertions | 5 |

**Total:** 9 API calls · ~2 min end-to-end
""")

    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Window size", "1,048,576 bp")
    col2.metric("Chromosome", "chr10")
    col3.metric("Cell type", "Mouse ESC")
    col4.metric("Total API calls", "9")

    if os.path.exists(MEDIA["exp1"]):
        st.image(MEDIA["exp1"], caption="Exp 1 — WT baseline", use_container_width=True)


elif page == "Exp 1 — Baseline":
    st.title("Experiment 1 — Baseline Validation")
    st.markdown("""
Confirms AlphaGenome correctly predicts the Unc5b CTCF cluster.
- **Panel 1:** Predicted contact heatmap
- **Panel 2:** CTCF ChIP-seq signal with peak annotations
- **Panel 3:** Insulation score track
""")

    if os.path.exists(MEDIA["exp1"]):
        st.image(MEDIA["exp1"], use_container_width=True)
    else:
        st.warning("Figure not yet generated. Run `make run-synthesis-exp1` or use the **Run Experiments** page.")

    if cache_exists:
        st.subheader("WT Cache Inspector")
        wt = load_wt_cache()
        if wt:
            col1, col2, col3 = st.columns(3)
            col1.metric("Contact map", f"{wt['n_bins']}×{wt['n_bins']}")
            col2.metric("Resolution", f"{wt['resolution']:,.0f} bp/bin")
            col3.metric("CTCF signal bins", str(len(wt["ctcf_signal"])))

            st.subheader("CTCF Signal")
            xbins = np.linspace(WIN_START / 1e6, WIN_END / 1e6, len(wt["ctcf_signal"]))
            fig, ax = plt.subplots(figsize=(10, 2.5))
            ax.fill_between(xbins, wt["ctcf_signal"], alpha=0.6, color="steelblue")
            ax.set_xlabel("Genomic position (Mb)")
            ax.set_ylabel("CTCF signal")
            ax.set_title("CTCF ChIP-seq signal — Unc5b window")
            st.pyplot(fig)
            plt.close(fig)

            st.subheader("Insulation Score")
            ins = compute_insulation(wt["contact_mat"], wt["resolution"])
            xbins2 = np.linspace(WIN_START / 1e6, WIN_END / 1e6, wt["n_bins"])
            fig2, ax2 = plt.subplots(figsize=(10, 2.5))
            ax2.plot(xbins2, ins, color="darkgreen", lw=1.5)
            ax2.set_xlabel("Genomic position (Mb)")
            ax2.set_ylabel("Insulation score")
            st.pyplot(fig2)
            plt.close(fig2)


elif page == "Exp 2 — Deletions":
    st.title("Experiment 2 — CTCF Deletion Series")
    st.markdown("""
Three loss-of-function variants probe the dose-dependence of boundary strength:

| Variant | Description |
|---------|-------------|
| **del_2ctcf** | Delete 2 weakest CTCF peaks |
| **del_4ctcf** | Delete all 4 weakest peaks |
| **flip_orient** | Flip strongest peak to divergent orientation |
""")

    tab1, tab2 = st.tabs(["Contact Panels", "Δ Insulation"])
    with tab1:
        if os.path.exists(MEDIA["exp2_panels"]):
            st.image(MEDIA["exp2_panels"], use_container_width=True)
        else:
            st.warning("Run `make run-synthesis-exp2` to generate this figure.")
    with tab2:
        if os.path.exists(MEDIA["exp2_delta"]):
            st.image(MEDIA["exp2_delta"], use_container_width=True)
        else:
            st.warning("Run `make run-synthesis-exp2` to generate this figure.")


elif page == "Exp 3 — Synthetic Insulator":
    st.title("Experiment 3 — Synthetic Stronger Insulator")
    st.markdown("""
Convergent CTCF pairs are inserted into the lowest-signal gap in the cluster.
Five designs are ranked by Δ insulation (positive = stronger than WT).

**CTCF motifs used:**
""")
    col1, col2 = st.columns(2)
    col1.code(f"FWD: {CTCF_FWD}")
    col2.code(f"REV: {CTCF_REV}")

    st.markdown("**Design preview:**")
    for n in [1, 2, 3]:
        seq = build_ctcf_pair_seq(n)
        st.code(f"{n} pair(s) ({len(seq)} bp): {seq[:60]}{'...' if len(seq) > 60 else ''}")

    if os.path.exists(MEDIA["exp3"]):
        st.image(MEDIA["exp3"], use_container_width=True)
    else:
        st.warning("Run `make run-synthesis-exp3` to generate this figure.")


elif page == "Run Experiments":
    st.title("Run Experiments")
    st.info(
        "These buttons trigger live AlphaGenome API calls. "
        "Ensure `ALPHA_GENOME_API_KEY` is set in your `.env` file."
    )

    st.subheader("Experiment 1 — Baseline")
    if st.button("▶ Run Exp 1 (1 API call, ~15 s)"):
        with st.spinner("Running exp1_baseline.py ..."):
            sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))
            try:
                import importlib, synthesis.exp1_baseline as e1
                importlib.reload(e1)
                model, dna_client = e1.load_model()
                wt_data = e1.predict_wildtype(model, dna_client)
                peaks   = e1.find_ctcf_peaks(wt_data["ctcf_signal"])
                e1.save_cache(wt_data)
                e1.plot_baseline(wt_data, peaks)
                st.success(f"Done! Figure saved to {e1.OUT_PNG}")
                st.image(e1.OUT_PNG, use_container_width=True)
            except Exception as exc:
                st.error(f"Failed: {exc}")

    st.subheader("Experiment 2 — CTCF Deletions")
    if st.button("▶ Run Exp 2 (3 API calls, ~45 s)"):
        with st.spinner("Running exp2_ctcf_deletions.py ..."):
            try:
                import importlib, synthesis.exp2_ctcf_deletions as e2
                importlib.reload(e2)
                e2.main()
                st.success("Done!")
                for path in [e2.OUT_PANELS, e2.OUT_DELTA]:
                    if os.path.exists(path):
                        st.image(path, use_container_width=True)
            except Exception as exc:
                st.error(f"Failed: {exc}")

    st.subheader("Experiment 3 — Synthetic Insulator")
    if st.button("▶ Run Exp 3 (5 API calls, ~75 s)"):
        with st.spinner("Running exp3_synthetic_insulator.py ..."):
            try:
                import importlib, synthesis.exp3_synthetic_insulator as e3
                importlib.reload(e3)
                e3.main()
                st.success("Done!")
                if os.path.exists(e3.OUT_PNG):
                    st.image(e3.OUT_PNG, use_container_width=True)
            except Exception as exc:
                st.error(f"Failed: {exc}")

    st.subheader("Generate Report")
    if st.button("📄 Generate synthesis_report.html"):
        try:
            import importlib, synthesis.synthesis_report as sr
            importlib.reload(sr)
            sr.main()
            st.success("synthesis_report.html generated!")
            with open("synthesis_report.html") as f:
                st.download_button(
                    "⬇ Download synthesis_report.html",
                    f.read(), file_name="synthesis_report.html", mime="text/html",
                )
        except Exception as exc:
            st.error(f"Failed: {exc}")

# Work Order 04: Visual Interpretation & Causal Dashboard

**Status:** Proposed (Priority: High)
**Owner:** HG-DT Core Implementation Team
**Context:** Interpretation & Visual Insight (Phase 4)

## 1. Objective
Build an interactive, modular Streamlit dashboard that visualizes the multi-scale consequences of DNA modifications—from 1D browser tracks to 2D chromatin maps and 3D protein structures—providing "Mechanistic Attribution" for all predictions.

## 2. Biological Background & Context
Visualization is the final step in turning a "Digital Twin" predicted into actionable biology. To understand a causal shift, the user needs to see the change at every level of the genome's hierarchy.

### Key Biological Concepts:
-   **1D Genomic Tracks:** Linear browser views for accessibility (ATAC), TF binding (ChIP), and expression (RNA-seq). Critical for spotting element-level changes.
-   **2D Chromatin Maps:** Hi-C triangular heatmaps representing the 3D folding of DNA. Comparing Ref vs. Mut heatmaps shows loop/TAD reorganization.
-   **3D Protein Structure:** The final physical product. Structural changes indicate functional loss or gain.
-   **Causal Attribution:** Linking a DNA mod (e.g., enhancer deletion) to a distal effect (e.g., Gene Y downregulation). Requires understanding enhancer–promoter loops.

### Where to Find More Information:
-   **pyGenomeTracks:** [pyGenomeTracks Documentation](https://pygenometracks.readthedocs.io/en/latest/) for 1D browser plotting.
-   **HiGlass:** [HiGlass: Web-based Visual Exploration](https://higlass.io/) for high-resolution interactive contact maps.
-   **NGLview:** [NGLView Widget Documentation](https://nglviewer.org/nglview/latest/) for 3D molecular visualization.
-   **Streamlit:** [Streamlit Main Documentation](https://docs.streamlit.io/) for fast Python-backed UI building.

## 3. Implementation Requirements

### A. Modular Viz Library
Build `src/hg_dt/viz/` containing:
1.  **`tracks_plotter.py`:** Wrapper for `pyGenomeTracks` to generate 1D "Reference vs. Mutant" browser plots (stored as PNG or interactive SVG).
2.  **`hic_plotter.py`:** Triangular Hi-C heatmap generator for the Ref vs. Mut 3D folding, using `cooler` and `cooltools`.
3.  **`protein_viz.py`:** Streamlit component for `NGLview` / `py3Dmol` to rotate and inspect 3D protein folding.

### B. Unified Streamlit Application
Develop `app.py` at the root with these tabs:
-   **Specification:** Input the DNA modification (chrom, locus, edit).
-   **Genome Tracks:** Show 1D deltas (accessibility/expression).
-   **3D Organization:** Show triangular contact map deltas + loop reorganization.
-   **Protein Structure:** View 3D folding (Ref vs. Mut).
-   **Trajectory Animation:** If the simulator is active, animate the differentiation time-course.
-   **Mechanistic Attribution:** A central text box using the "Causal Wrapper" to explain the phenotype.

### C. Automated Attribution Engine
Build `src/hg_dt/analyze/attribution.py` to:
-   Generate a "Mechanistic Insight" text string summarizing the multi-scale delta.
-   Example: *"This deletion removes ccRE EH38E1800647 → weakens enhancer–promoter loop → accessibility ↓28% → expression ↓35%."*

## 4. Implementation Steps
1.  **Develop linear plotter:** Integrate `pyGenomeTracks`.
2.  **Develop contact map plotter:** Build the triangular heatmap delta visualizer.
3.  **Integrate 3D Protein Viewer:** Configure `py3Dmol` in Streamlit.
4.  **Execute Phase 4 Dashboard:** Wire the `app.py` UI to all analytical modules.
5.  **Test Case:** Full sequence-to-structure-to-viz for the *TAL1* deletion and verify the dashboard successfully highlights the distal impact.

## 5. Success Criteria
-   Dashboard launches with single command `streamlit run app.py`.
-   Users can interactively compare Ref vs. Mut structures and contact maps side-by-side.
-   "Mechanistic Summary" correctly identifies the disrupted enhancer loop in the *TAL1* experiment.

---
*End of Framework Implementation Work Orders.*

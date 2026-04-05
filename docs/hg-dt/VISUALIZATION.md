# Visualization & Interpretation Layer (v1.1)

The **Visualization & Interpretation Layer (Step 10)** transforms raw model predictions into biologically meaningful insights and publication-ready figures.

## 1. Unified Dashboard (Streamlit MVP)

A central, interactive dashboard integrates all analytical results:

-   **Selector Tabs:** Genome Tracks | Contact Maps | 3D Trajectory | Protein 3D | Summary.
-   **Interactivity:** Sliders to adjust differentiation time-points or switch between K562 and GM12878 cell lines.
-   **Comparison Mode:** Side-by-side or overlaid "Reference vs. Mutant" visualizations.

## 2. Linear Genomic Tracks (1D Visualization)

-   **Tools:** `pyGenomeTracks`.
-   **Outputs:** Tracks for chromatin accessibility (ATAC/DNase), H3K27ac marks, TF binding, and predicted gene expression (RNA-seq).
-   **Interpretation:** Identifies the precise regulatory element disrupted by a DNA edit (e.g., loss of a promoter-proximal enhancer signal).

## 3. Chromatin Organization (2D Heatmaps)

-   **Tools:** `cooler`, `cooltools`, `HiCExplorer`.
-   **Outputs:** Triangular Hi-C heatmaps, TAD boundary calls, and differential contact maps.
-   **Predicted Loops:** Explicit loops overlaid as arcs, highlighting predicted loop strength deltas (↑/↓).
-   **3D Trajectory Animation:** An animated GIF showing TAD reorganization during differentiation.

## 4. Protein Structure & Impact (3D Molecular View)

-   **Tools:** `NGLview`, `py3Dmol`.
-   **Outputs:** 3D protein structures colored by pLDDT confidence (Blue: High, Red: Low).
-   **PAE Matrices:** Visualizes the Predicted Aligned Error (PAE) to identify domain shifts or lost interfaces.
-   **Functional Annotation:** Highlights catalytic sites or interaction domains mapped from UniProt/PDB.

## 5. Automated Mechanistic Attribution

Every DNA edit generates a text-based "Mechanistic Insight" summary:

> **Example:** "Deletion at `chrX:123,456` removes an ENCODE-mapped enhancer (EH38E1800647). AlphaGenome predicts a **42% reduction in contact probability** with the target *Gene Y* promoter and a subsequent **35% drop in total transcript abundance**. The resulting protein retains its primary fold but loses the C-terminal interaction domain (pLDDT shift from 88 to 45)."

## 6. Publication-Ready Exports

-   **Figures:** Vectorized PDFs or high-resolution TIFF files.
-   **Coordinates:** Export predicted tracks to BigWig or BED formats for use in third-party browsers (IGV / UCSC).
-   **Structures:** Export PDB/mmCIF files with predicted metadata.

---
*Reference: Visualization strategies for the Human Genome Digital Twin Prototype (April 2026).*

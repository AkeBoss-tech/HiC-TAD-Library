# HG-DT Transition Plan: From Mouse TADs to Human Digital Twin

This plan outlines the specific steps to transition the current `HiC-TAD-Library` (focused on mouse Micro-C analysis) to the comprehensive **Human Genome Digital Twin (HG-DT)** framework.

## Current Project Assessment
- **Status:** Functional mouse Hi-C analysis library with 3D polymer simulation and TAD boundary detection.
- **Data:** Primarily `mouse_microc.mcool` (mm10).
- **Architecture:** Procedural visualization scripts calling core analytical modules in `src/`.
- **Gap:** Lacks hg38 context, integrated DNA->RNA->Protein causal pipeline, and automated human regulatory element mapping (ENCODE cCREs).

## Transition Phases

### Phase 1: Human (hg38) Infrastructure & Annotation (Weeks 1-2)
*Transitioning project context from mm10 to hg38.*
1.  **HG-DT Package Setup:** Initialize `src/hg_dt` module to separate new causal logic from existing TAD utilities.
2.  **Annotation Layer:** Build the `ReferenceContextBuilder` to ingest GENCODE and ENCODE datasets.
3.  **Data Acquisition:** Script the download of hg38 FASTA, GENCODE v49, and ENCODE cCREs.

### Phase 2: System 2 - AlphaGenome & Regulatory 3D (Weeks 3-4)
*Bridging existing Hi-C logic with mutation-aware 3D prediction.*
1.  **AlphaGenome Integration:** Implement the `AlphaGenomeConnector` for dual-pass (Ref/Mut) 1 Mb predictions.
2.  **3D Organization Delta:** Refactor `src/tad_boundaries.py` to support comparing predicted Ref vs. Mut insulation and TAD boundaries.
3.  **Accessibility Scanning:** Implement accessibility track analysis to identify regulator silencing.

### Phase 3: System 1 - The Central Dogma Pipeline (Weeks 5-7)
*Implementing the core simulation from sequence to structure.*
1.  **Splicing & RNA Prediction:** Use sequence models to forecast mRNA isoforms from DNA edits.
2.  **Translation Engine:** Implementation of rule-based Biopython translation logic.
3.  **Protein Structure Profiler:** Chaining transcripts to AlphaFold 3 (BioNeMo/Local).

### Phase 4: Visualization, Causal UI & Evaluation (Weeks 8-12)
*Finalizing the research platform.*
1.  **Interpretable Dashboards:** Build the Streamlit dashboard for multi-scale visualization.
2.  **Differentiation Simulator:** Add the dynamic contact map animation layer.
3.  **Benchmark Suite:** Evaluate on known human variants (eQTL/ClinVar).

---
*Work orders in `spec/work_orders/` detail the exact implementation for each step.*

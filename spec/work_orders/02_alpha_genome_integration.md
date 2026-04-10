# Work Order 02: AlphaGenome Integration & 3D Organization Delta

**Status:** Proposed (Priority: High)
**Owner:** HG-DT Core Implementation Team
**Context:** System 2: Non-Coding Regulatory (Phase 1)

## 1. Objective
Enable the HG-DT platform to predict the regulatory and 3D organization consequences of DNA modifications using **AlphaGenome (Google DeepMind)** and **Evo 2 (Arc Institute)** models.

## 2. Biological Background & Context
This system models the "Non-Coding Regulatory" landscape. Mutations in enhancers or promoters can silence genes without changing the protein sequence directly. Understanding 3D organization is the key to unlocking the mechanism.

### Key Biological Concepts:
-   **Chromatin Accessibility:** Regions of the genome where the DNA is "open" and accessible to transcription factors (TFs). Measured by ATAC-seq or DNase-seq. Silencing an enhancer often reduces accessibility.
-   **Hi-C Contact Maps (3D Folding):** Heatmaps showing how frequently different parts of the genome physically touch. High-frequency loops often connect enhancers to distal promoters.
-   **TAD Boundaries:** Large structural domains (Topologically Associating Domains) that insulate regulation. A mutation that disrupts a Boundary can lead to misregulation.
-   **Variant Effects:** The magnitude of change (delta) in accessibility or contacts between a wild-type and mutated sequence.

### Where to Find More Information:
-   **AlphaGenome:** [AlphaGenome Publication (Nature 2026)](https://www.nature.com/articles/s41586-026-00000-0)
-   **Evo 2:** [Evo 2 Research (Arc Institute 2026)](https://arcinstitute.org/tools/evo)
-   **Hi-C/TADs:** See [Dixon et al. (2012)](https://www.nature.com/articles/nature11082) for TAD discovery and [Rao et al. (2014)](https://www.cell.com/cell/fulltext/S0092-8674(14)01497-4) for high-resolution loop mapping.
-   **Variant Effects:** [Enformer (2021)](https://www.nature.com/articles/s41592-021-01252-x) provides the foundation for sequence-to-track prediction.

## 3. Implementation Requirements

### A. AlphaGenomeConnector Class
Implement `src/hg_dt/models/alphagenome.py` to:
1.  **Interface with AlphaGenome API:** Call the model with a 1 Mb sequence (Ref/Mut) and a list of output tracks.
2.  **Predict All Relevant Tracks:** Accessibility (ATAC), TF binding, Histone marks, CAGE (transcription initiation), and Hi-C maps.
3.  **Handle Multi-Resolution:** Extract contact maps at the native ~2 kb resolution (or interpolate to match existing library utilities).

### B. Track Delta Engine
Build `src/hg_dt/analyze/deltas.py` to:
-   **Normalize Predictions:** Scale Ref/Mut outputs to ensure comparable units.
-   **Compute Quantitative Deltas:**
    -   `accessibility_delta`: $\log_2(Ref/Mut)$ across the 1 Mb window.
    -   `contact_delta`: Difference heatmap for predicted loops.
    -   `expression_delta`: Predicted RNA abundance change for target promoters.
-   **Identify "Silenced" Elements:** Find ENCODE cCREs where predicted accessibility drops by >50%.

### C. Refactoring TAD Logic
-   Update `src/tad_boundaries.py` to support `compute_insulation_delta(ref_map, mut_map)` and `compare_tad_boundaries(ref_tads, mut_tads)`.
-   Enable the 3D polymer simulation (`src/polymer_sim.py`) to run on predicted Ref vs. Mut contact maps.

## 4. Implementation Steps
1.  **Setup API Connectivity:** Configure `ALPHA_GENOME_API_KEY` in `.env`.
2.  **Develop `AlphaGenomeConnector`:** Build the wrapper for the DeepMind client.
3.  **Implement `analyze.compute_deltas`:** Create the logic for track-wise comparison.
4.  **Develop Loop Identification:** Implement `hg_dt.analyze.find_distal_loops(contact_map)` to identify promoter-enhancer pairs.
5.  **Integration Test:** Predict the 3D organizational shift for a known *MYC* enhancer deletion and compare against literature results.

## 5. Success Criteria
-   `AlphaGenomeConnector` consistently returns all regulatory tracks for a 1 Mb sequence.
-   `compute_track_deltas` correctly identifies significant loop dissolution in the *MYC* locus after mock enhancer deletion.
-   System successfully outputs a quantitative "Likelihood of Expression Loss" for target genes.

---
*Next Task: Work Order 03 - The DNA-to-Protein Pipeline.*

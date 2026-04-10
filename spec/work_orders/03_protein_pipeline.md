# Work Order 03: The DNA-to-Protein Pipeline

**Status:** Proposed (Priority: High)
**Owner:** HG-DT Core Implementation Team
**Context:** System 1: Gene/Protein Impact (Phase 2)

## 1. Objective
Enable the HG-DT platform to go from a DNA modification (SNP, deletion, insertion) to a functional protein predictively. This covers transcription, alternative splicing, rule-based translation, and 3D protein structure.

## 2. Biological Background & Context
The "Central Dogma" of biology: DNA is transcribed into mRNA, and mRNA is translated into protein. A mutation can affect any stage.

### Key Biological Concepts:
-   **Transcription:** RNA Polymerase reads DNA to create a single-stranded pre-mRNA.
-   **Alternative Splicing:** The spliceosome removes introns and joins exons. Different combinations of exons (isoforms) can create different proteins from the same gene.
-   **Translation:** Ribosomes read mRNA as codons (3 bases) to link amino acids. A single-base deletion can cause a **frameshift**, completely changing the protein sequence.
-   **Protein Folding (3D):** The amino acid chain folds into a 3D shape. A single mutation can destabilize a fold or break an interaction domain.
-   **pLDDT (Predicted Local Distance Difference Test):** A per-residue confidence metric (0-100) from AlphaFold. 90+ is very high; <50 is unreliably folded.

### Where to Find More Information:
-   **Central Dogma:** See [Alberts et al., *Molecular Biology of the Cell*](https://www.ncbi.nlm.nih.gov/books/NBK21054/) or [Central Dogma Overview (Khan Academy)](https://www.khanacademy.org/science/biology/gene-expression-central-dogma/central-dogma-transcription-translation/a/the-central-dogma-of-molecular-biology).
-   **AlphaFold 3:** [AlphaFold 3 Publication (Nature 2024/2026)](https://www.google.com/search?q=AlphaFold+3+DeepMind+Nature) and [AlphaFold Protein Structure Database](https://alphafold.ebi.ac.uk/).
-   **Splicing Prediction:** See [SpliceAI (2019)](https://www.cell.com/cell/fulltext/S0092-8674(18)31629-5) for sequence-to-splice-junction prediction foundation.
-   **Proteomics:** Consult [UniProt](https://www.uniprot.org/) for thousands of verified human protein structures and functional annotations.

## 3. Implementation Requirements

### A. Transcript/mRNA Pipeline
Implement `src/hg_dt/translate/transcript.py` to:
1.  **Extract Exons:** Use GENCODE annotations from `ReferenceContextBuilder`.
2.  **Predict Split Effects:** Use sequence models (Evo 2 or splicing tracks from AlphaGenome) to forecast if a mutation causes exon skipping or truncation.
3.  **Assemble mRNA Isoforms:** Construct the mature mRNA sequence string.

### B. Translation Engine
Implement `src/hg_dt/translate/translator.py` using `Biopython`:
-   **Identify Open Reading Frames (ORFs):** Scan mRNA for the start codon (**AUG**) and stop codons (UAG, UAA, UGA).
-   **Codon-to-AA Translation:** Map the genetic code to an amino acid string.
-   **Detect Frameshifts:** Flag if a DNA modification shifts the triplet reading frame.

### C. Protein Structure Profiler
Implement `src/hg_dt/models/protein.py` to:
-   **Call AlphaFold 3 / ESMFold:** Send the Ref/Mut amino acid sequences to the folding engine (via BioNeMo or local instance).
-   **Structural Impact Scanning:** Compute structural deltas between wild-type and mutant (e.g., lost domains or destabilized pockets).

## 4. Implementation Steps
1.  **Develop `transcript.predict_isoforms()`:** Predict the mature mRNA(s) following DNA edits.
2.  **Develop `translator.translate()`:** Create the logic for rule-based codon translation.
3.  **Orchestrate Protein Folding:** Integrate BioNeMo’s AlphaFold 3 NIM/API for Ref vs. Mut structure prediction.
4.  **Develop `analyze.compute_structure_impact()`:** Implement metrics to identify fold collapse (e.g., pLDDT drops).
5.  **Integration Test:** Predict the protein structure shift for a known *BRCA1* pathogenic SNP and compare against ClinVar and AlphaFold DB.

## 5. Success Criteria
-   `translate_to_protein` correctly handles frameshifts for 1 bp DNA deletions.
-   Successfully predicts structural collapse (pLDDT drop) for known truncating variants.
-   Platform can output a 3D `.pdb` file for both Ref and Mut proteins.

---
*Next Task: Work Order 04 - Visual Interpretation & Causal Dashboard.*

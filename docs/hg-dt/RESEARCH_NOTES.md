# Biological Context & Data Sources

This document provides foundational research on the **central dogma of molecular biology** (DNA → RNA → Protein) and details the critical **data sources** for the HG-DT project (Human Genome Digital Twin).

## 1. The Multi-Step Pipeline: DNA → RNA → Protein

A "full simulated cell" pipeline must model each stage of the central dogma to accurately predict the consequences of a DNA edit.

### A. DNA → RNA (Transcription & Processing)
-   **Initiation:** RNA polymerase II binds to a promoter, often aided by enhancers.
-   **Processing:**
    -   **5' Capping:** Modified guanine for stability.
    -   **Splicing:** Removal of introns; alternative splicing combines exons into multiple isoforms.
    -   **3' Poly-A Tail:** Stability and export signal.
-   **Modeling with ML:**
    -   **AlphaGenome (2026):** Predicts gene expression levels, transcription initiation, and splice site usage from 1 Mb of DNA.
    -   **Borzoi (2025):** Predicts cell-type-specific RNA-seq coverage.
    -   **Evo 2 (2026):** Predicts mutational impacts on alternative splicing and RNA secondary structure.

### B. RNA → Protein (Translation & Folding)
-   **Translation:** Ribosomes read mRNA codons (sets of 3 nucleotides) to link amino acids.
-   **Efficiency:** Influenced by 5' UTR sequence, codon usage, and mRNA secondary structure.
-   **Modeling with ML:**
    -   **Biopython:** Rule-based sequence translation.
    -   **RiboNN (2025):** Predicts protein yield (translation efficiency) across cell types.
    -   **AlphaFold 3:** Predicts the 3D folded structure and interactions of the resulting protein.

## 2. Definitive Human Genomic Data Sources (hg38)

To map gene locations and regulatory elements, the project leverages the following standard repositories:

| Resource | Primary Purpose | Link/Format |
| :--- | :--- | :--- |
| **GENCODE Release 49** | Comprehensive gene annotations (Gold Standard). | [Gencode Human](https://www.gencodegenes.org/human/) (GTF) |
| **ENCODE Project** | Registry of candidate Cis-Regulatory Elements (cCREs). | [Screen Registry](https://screen.wenglab.org/downloads) (BED) |
| **UCSC Genome Browser** | Visualization context and bulk track downloads (H3K27ac, etc.). | [hg38 Gateway](https://genome.ucsc.edu/cgi-bin/hgGateway) |
| **4D Nucleome (4DN)** | Uniformly processed Hi-C and Micro-C contact maps. | [4DN Data Portal](https://data.4dnucleome.org/) |

### Key Regulatory Concepts
-   **cCREs:** Candidate Cis-Regulatory Elements, mapping millions of enhancers, promoters, and TF binding sites across 1,888 cell types.
-   **Hi-C Interactions:** Essential for identifying how distant enhancers (mapped by ENCODE) loop to meet target promoters.

## 3. Practical Data Preparation (One-Liners)

Use these commands to download the baseline annotations for System 1 and System 2:

```bash
# GENCODE GTF (Genes, Transcripts, Exons)
wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

# ENCODE cCREs (Promoters, Enhancers, TF motifs)
wget https://screen.wenglab.org/downloads/GRCh38-cCREs.bed.gz
```

---
*Reference Context: Central Dogma and Genomic Repositories (HG-DT v0.1 Research).*

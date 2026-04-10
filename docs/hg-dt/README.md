# Human Genome Digital Twin (HG-DT) Prototype

Welcome to the **Human Genome Digital Twin (HG-DT)** project documentation. This initiative at the Kwan Lab aims to build a research-grade, causal multi-omics platform to simulate the effects of DNA modifications on cellular states.

## Project Overview

The HG-DT platform is designed to answer critical questions about genomic intervention: **"If I edit this specific section of DNA, what happens to the 3D chromatin organization, gene expression, and the resulting proteins?"**

Unlike descriptive models that merely correlate sequence with function, HG-DT focuses on **interventional prediction**, providing a virtual laboratory for understanding the consequences of SNPs, indels, large deletions, and regulatory mutations.

## Core Systems

The prototype (v0.1) integrates two primary analytical systems:

1.  **Gene-Centric Protein Impact System:** Maps DNA sections to genes and regulators to predict how deletions or mutations alter transcript isoforms and downstream protein structure/function.
2.  **Regulatory 3D-Organization System:** Leverages AlphaGenome and Evo 2 to predict changes in chromatin accessibility and 3D contact maps (Hi-C) following DNA edits, and correlates these structural shifts with expression levels.

## Documentation Structure

- [**SPECIFICATION.md**](./SPECIFICATION.md): Detailed system architecture, technical stack, research plan, and milestones.
- [**RESEARCH_NOTES.md**](./RESEARCH_NOTES.md): Biological foundations (DNA → RNA → Protein) and comprehensive data sources (ENCODE, GENCODE, etc.).
- [**VISUALIZATION.md**](./VISUALIZATION.md): Strategies for interpreting and visualizing the system's products, from contact maps to 3D protein structures.

## Key Technologies

- **DNA Foundation Models:** AlphaGenome, Evo 2, Borzoi.
- **Protein Structure:** AlphaFold 3, ESMFold.
- **Orchestration:** NVIDIA BioNeMo, Python.
- **Bioinformatics:** pyGenomeTracks, cooler, GENCODE, ENCODE Screen Registry.

---
*Developed for the Kwan Lab 3D Genome Research Program.*

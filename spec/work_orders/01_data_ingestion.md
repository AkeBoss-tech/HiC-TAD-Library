# Work Order 01: hg38 Data Infrastructure & ReferenceContextBuilder

**Status:** Proposed (Priority: High)
**Owner:** HG-DT Core Implementation Team
**Context:** Environment & Core Data Pipelining (Phase 0)

## 1. Objective
Establish the **Human hg38 (GRCh38)** computational environment and build a centralized `ReferenceContextBuilder` that aggregates genome sequence, gene annotations (GENCODE), and regulatory element registries (ENCODE).

## 2. Biological Background & Context
For a "Digital Twin" to be research-grade, it MUST be grounded in the latest human reference assembly and high-quality annotations.

### Key Biological Concepts:
-   **hg38 (GRCh38):** The standard reference assembly for the human genome. It provides the coordinate system for all genomic data.
-   **GENCODE v49:** The "Gold Standard" for human gene annotations. Unlike basic transcripts, GENCODE provides comprehensive mapping for exons, introns, isoforms, and non-coding RNAs.
-   **ENCODE Screen Registry (cCREs):** Maps million of Candidate Cis-Regulatory Elements (cCREs). Enhancers and promoters are often located kilobases away from their target genes, making their 3D context (mapped via Hi-C) critical.

### Where to Find More Information:
-   **DNA Sequence:** [UCSC Genome Browser hg38](https://genome.ucsc.edu/cgi-bin/hgGateway?db=hg38)
-   **Gene Models:** [GENCODE Main Site](https://www.gencodegenes.org/human/)
-   **Regulatory Elements:** [ENCODE Screen cCREs Documentation](https://screen.wenglab.org/docs)
-   **File Formats:** [GTF Specification](https://m.ensembl.org/info/website/upload/gff.html) (for GENCODE) and [BED Specification](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) (for ENCODE).
-   **Hi-C & TADs:** See [4D Nucleome Data Portal](https://data.4dnucleome.org/) for uniformly processed human Hi-C.

## 3. Implementation Requirements

### A. Data Ingestion Scripts
Create a `scripts/setup_hg38_env.sh` to download:
1.  **hg38 FASTA:** `rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr*.fa.gz .` or from `bigZips`.
2.  **GENCODE v49 GTF:** `wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz`
3.  **ENCODE cCREs BED:** `wget https://screen.wenglab.org/downloads/GRCh38-cCREs.bed.gz`

### B. ReferenceContextBuilder Class
Implement `src/hg_dt/data/builder.py` with these capabilities:
-   **`pyfaidx` Integration:** Fast, random-access sequence extraction from local `.fa` files.
-   **`pyranges` for Annotations:** Efficiently query the GENCODE GTF and ENCODE BED files to find all features overlapping a given 1 Mb window.
-   **Window Expansion:** Given an input coordinate (SNP or deletion), expand to a **1 Mb window** centered on the edit (matching AlphaGenome's input requirement).
-   **Mutant Sequence Generation:** Programmatic generation of the "edited" FASTA string by applying the modification (SNP replacement, deletion, or insertion) to the reference sequence.

## 4. Implementation Steps
1.  **Setup Directory:** `mkdir -p data/reference/hg38` and `mkdir -p src/hg_dt/data`.
2.  **Develop Data Downloader:** Create `scripts/fetch_reference_data.py` to automate the downloads above.
3.  **Develop `ReferenceContextBuilder`:**
    -   Parse GTF (GENCODE) to build a gene/transcript lookup table.
    -   Parse BED (ENCODE) to build a regulatory element lookup.
    -   Implement `get_context(chrom, start, end, edit_type, new_sequence)` which returns:
        -   `ref_seq`: 1 Mb reference FASTA string.
        -   `mut_seq`: 1 Mb mutant FASTA string.
        -   `annotations`: List of genes and cCREs in the window.
4.  **Test Case:** Verify results against a known 1 Mb region on chromosome 17 (e.g., surrounding *BRCA1*).

## 5. Success Criteria
-   `ReferenceContextBuilder.get_context()` executes in < 200ms.
-   Correctly identifies all exons and enhancers within the *BRCA1* 1 Mb neighborhood.
-   Accurately generates the mutant FASTA for a mock 5 kb enhancer deletion.

---
*Next Task: Work Order 02 - AlphaGenome & Regulatory 3D.*

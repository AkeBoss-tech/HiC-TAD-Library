# LinkPrep Integration Plan

This folder tracks the work needed to support Dovetail LinkPrep data cleanly in this repository.

## Goal

Take Dovetail LinkPrep paired-end FASTQs through:

1. preprocessing to `.pairs`, BAM, and `.mcool`
2. QC checks against Dovetail whole-genome recommendations
3. downstream analysis for:
   - compartments
   - TADs
   - loops
   - SV / CNV
   - phasing
4. a clearer separation between:
   - Dovetail-recommended external tools
   - this repo's internal exploratory analysis code

## Current State

What already exists in this repo:

- Raw-data preprocessing script:
  - [`scripts/preprocess_microc.sh`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/scripts/preprocess_microc.sh)
- Internal downstream analysis modules:
  - [`src/tad_boundaries.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/src/tad_boundaries.py)
  - [`src/compartments.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/src/compartments.py)
  - [`src/loops.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/src/loops.py)
  - [`src/sv_detection.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/src/sv_detection.py)
  - [`src/phasing.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/src/phasing.py)
- Visualization / demo entrypoints:
  - [`visualize_compartments.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/visualize_compartments.py)
  - [`visualize_loops.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/visualize_loops.py)
  - [`analyze_sv_cnv.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/analyze_sv_cnv.py)
  - [`analyze_phasing.py`](/Users/akashdubey/Documents/CodingProjects/HiC-TAD-Library/analyze_phasing.py)

What is missing right now on this machine:

- `bwa`
- `samtools`
- `pairtools`
- `pairix`
- `bgzip`
- LinkPrep FASTQ inputs
- `mm10.chromsizes`

What is still missing in the repo itself:

- a LinkPrep-specific runbook
- a cleaner config layer for preprocessing inputs
- automated QC summarization
- explicit integration of Dovetail-recommended downstream tools
- tests around the preprocessing path

## Dovetail Alignment

The preprocessing script already follows the documented Dovetail-style pattern:

1. `bwa mem -5SP -T0`
2. `pairtools parse`
3. `pairtools sort`
4. `pairtools dedup`
5. `pairtools split`
6. `pairix`
7. `cooler cload pairix`
8. `cooler zoomify --balance`

That is consistent with the Dovetail docs:

- [Data processing overview](https://dovetail-analysis.readthedocs.io/en/latest/)
- [FASTQ to final valid pairs bam](https://dovetail-analysis.readthedocs.io/en/latest/data_processing/fastq_to_bam.html)
- [Generating contact matrices](https://dovetail-analysis.readthedocs.io/en/latest/data_processing/contact_map.html)
- [Whole-genome QC](https://dovetail-analysis.readthedocs.io/en/latest/whole_genome/qc.html)
- [Whole-genome feature discovery](https://dovetail-analysis.readthedocs.io/en/latest/whole_genome/feature_discovery.html)
- [Whole-genome variant detection](https://dovetail-analysis.readthedocs.io/en/latest/whole_genome/genetic_variant_detection.html)

## Recommended Work Plan

### Phase 1: Make preprocessing runnable

- Install required command-line tools in the analysis environment.
- Add a chromosome sizes file for the chosen reference.
- Replace placeholder FASTQ paths with real LinkPrep inputs.
- Run preprocessing end to end and save outputs under `data/processed/`.

Suggested install command:

```bash
conda activate hic-analysis
conda install -c bioconda bwa samtools pairtools pairix cooler htslib
```

Suggested chromosome sizes command:

```bash
cut -f1,2 data/raw/mm10.fa.fai > data/raw/mm10.chromsizes
```

### Phase 2: Add QC and provenance

- Parse `pairtools dedup` stats into a machine-readable QC summary.
- Record:
  - total pairs
  - duplicate rate
  - no-dup pairs
  - cis pairs >= 1 kb
  - estimated suitability for TAD / loop calling
- Save a QC summary markdown or JSON report next to the generated `.mcool`.

Target thresholds from Dovetail docs:

- no-dup cis pairs >= 1 kb > 40% of no-dup pairs
- >125M no-dup pairs for strong TAD / loop analysis

### Phase 3: Clarify downstream modes

Separate downstream analysis into two categories:

- `internal exploratory`
  - current Python implementations in `src/`
- `dovetail-recommended external`
  - TADs: Arrowhead / Juicer ecosystem
  - loops: Mustache
  - compartments: FAN-C
  - SVs: hic_breakfinder
  - CNV: CNVkit
  - SNVs / indels: DeepVariant

This repo should be explicit about when a result is:

- a custom approximation
- a Dovetail-style preprocessing result
- a Dovetail-recommended downstream call

### Phase 4: Add a LinkPrep-specific entrypoint

Create a dedicated runner for LinkPrep, for example:

- `linkprep/run_linkprep.sh`
- `linkprep/qc_summary.py`
- `linkprep/README.md`

That runner should:

1. validate required tools
2. validate required inputs
3. run preprocessing
4. generate QC summary
5. print exact output file locations

## Immediate TODO

- [ ] Install missing preprocessing tools.
- [ ] Add `data/raw/mm10.chromsizes`.
- [ ] Decide on real LinkPrep FASTQ locations.
- [ ] Run `scripts/preprocess_microc.sh` on a real sample.
- [ ] Capture QC output into a persistent report.
- [ ] Add a LinkPrep-specific wrapper script under this folder.
- [ ] Add documentation for internal-vs-recommended downstream analyses.

## Proposed Next Files

These are the next files worth adding:

- `linkprep/run_linkprep.sh`
- `linkprep/qc_summary.py`
- `linkprep/config.example.sh`

## Notes From Initial Audit

- The repo can analyze existing `.mcool` files today.
- The raw LinkPrep preprocessing path is scaffolded but not currently runnable in this environment because required tools are missing.
- The current downstream Python modules are useful for exploration, but they are not drop-in replacements for all Dovetail-recommended external tools.


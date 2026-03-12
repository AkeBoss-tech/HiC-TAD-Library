# LinkPrep Analysis

End-to-end analysis pipeline for **Dovetail LinkPrep™** proximity-ligation data.

## What LinkPrep measures

LinkPrep uses **Tn5 tagmentation** (uniform, sequence-agnostic fragmentation) followed by proximity ligation. This produces linked reads that encode:

- 3D genome structure (TADs, compartments, loops)
- Structural variants (translocations, inversions, deletions)
- Copy-number variants (from uniform coverage depth)
- Haplotype phase blocks (from linked heterozygous SNPs)

---

## Folder layout

```
linkprep/
├── config.py                        ← DATA PATHS & PARAMETERS — edit this for real data
├── run_analysis.py                  ← Master entry point
│
├── sample_data/
│   └── generate_sample_data.py      ← Creates synthetic .cool, .pairs, .vcf, coverage.bed
│
├── pipeline/
│   ├── preprocess.sh                ← FASTQ → .cool pipeline (real data)
│   └── qc_metrics.py                ← Pairs QC: valid-pair rate, cis/trans, P(s) curve
│
├── analysis/
│   ├── contact_analysis.py          ← TADs, compartments, loops from .cool
│   ├── sv_detection.py              ← Translocation, inversion, deletion calling
│   ├── cnv_analysis.py              ← Per-bin coverage → CNV segments
│   └── phasing.py                   ← VCF SNPs → haplotype blocks
│
├── visualization/
│   ├── plot_contacts.py             ← Hi-C heatmap, TAD overlay, A/B compartments
│   ├── plot_sv.py                   ← Translocation heatmap, inversion map
│   ├── plot_cnv.py                  ← Genome-wide & per-chrom CNV
│   └── plot_phasing.py              ← SNP density, block map, length distribution
│
└── figures/                         ← All output PNGs (auto-created)
```

---

## Quick start (sample data)

```bash
cd HiC-TAD-Library
conda activate hic-analysis

# Run everything — generates synthetic data automatically
python linkprep/run_analysis.py
```

Figures are written to `linkprep/figures/`.

---

## Switching to real data

1. Run the preprocessing pipeline on your FASTQs:

```bash
bash linkprep/pipeline/preprocess.sh \
    /path/to/reference.fa \
    /path/to/sample_R1.fastq.gz \
    /path/to/sample_R2.fastq.gz \
    my_sample \
    16      # threads
```

2. Update `linkprep/config.py`:

```python
COOL_FILE    = "/path/to/my_sample.cool"
PAIRS_FILE   = "/path/to/my_sample.valid.pairs.gz"
VCF_FILE     = "/path/to/my_sample.vcf"
COVERAGE_BED = "/path/to/my_sample_coverage.bed"
GENOME       = "hg38"   # or mm10

REGIONS = {
    "MyGene_Chr1": "chr1:100,000,000-105,000,000",
}
```

3. Re-run the analysis:

```bash
python linkprep/run_analysis.py --no-generate
```

---

## Running individual steps

```bash
# Only QC + contact analysis
python linkprep/run_analysis.py --steps qc contacts

# Only SV + CNV
python linkprep/run_analysis.py --steps sv cnv --no-generate

# Individual scripts
python linkprep/pipeline/qc_metrics.py
python linkprep/visualization/plot_contacts.py
python linkprep/visualization/plot_sv.py
python linkprep/visualization/plot_cnv.py
python linkprep/visualization/plot_phasing.py
```

---

## Output figures

| File | Contents |
|------|----------|
| `qc_distance_distribution.png` | Cis-pair P(s) curve and raw distance histogram |
| `qc_pair_types.png` | Trans / short-range cis / long-range cis breakdown |
| `contact_heatmap_<region>.png` | Hi-C heatmap + insulation track |
| `tads_<region>.png` | Contact heatmap with TAD boundary boxes |
| `compartments_<region>.png` | Heatmap + E1 eigenvector (A=red, B=blue) |
| `sv_translocation_heatmap.png` | Inter-chromosomal pair density |
| `sv_summary.png` | SV count by type |
| `sv_inversions_<chrom>.png` | Putative inversion spans |
| `cnv_genome_wide.png` | Genome-wide log2(obs/exp) with CBS segments |
| `cnv_segments_<chrom>.png` | Per-chromosome segmentation detail |
| `phasing_snp_density.png` | Heterozygous SNP density |
| `phasing_haplotype_blocks.png` | Phase-block map per chromosome |
| `phasing_block_lengths.png` | Block length distribution + N50 |

---

## Key QC thresholds (LinkPrep)

| Metric | Pass threshold |
|--------|---------------|
| Valid (UU) pair rate | ≥ 50% |
| Cis pair fraction | ≥ 60% |
| Long-range cis (≥1 kb) | ≥ 40% |
| Complexity (unique fraction) | ≥ 70% |

---

## Dependencies

All are already in the `hic-analysis` conda environment:

```
cooler  numpy  pandas  scipy  matplotlib
```

For the preprocessing shell script:
```
bwa  samtools  pairtools
```

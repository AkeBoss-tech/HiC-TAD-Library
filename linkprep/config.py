"""
Central configuration for LinkPrep analysis.

When you receive real data, update the DATA section below and re-run
`python run_analysis.py`.  Everything else adapts automatically.
"""

import os

# ---------------------------------------------------------------------------
# DATA SOURCES — update these when real data arrives
# ---------------------------------------------------------------------------

# Directory where all input files live
DATA_DIR = os.path.join(os.path.dirname(__file__), "sample_data", "output")

# Contact matrix (.cool or .mcool).  For .mcool append ::resolutions/<res>
COOL_FILE = os.path.join(DATA_DIR, "linkprep_sample.cool")

# Raw pairs file produced by pairtools (gzipped or plain)
PAIRS_FILE = os.path.join(DATA_DIR, "linkprep_sample.pairs.gz")

# Variant file (VCF or VCF.gz) — used for phasing
VCF_FILE = os.path.join(DATA_DIR, "linkprep_sample.vcf")

# Per-bin coverage BED (chrom  start  end  coverage)
COVERAGE_BED = os.path.join(DATA_DIR, "linkprep_sample_coverage.bed")

# Reference genome name (used in figure labels only)
GENOME = "mm10"

# ---------------------------------------------------------------------------
# ANALYSIS PARAMETERS
# ---------------------------------------------------------------------------

# Bin size used in the .cool file (bp)
BIN_SIZE = 10_000          # 10 kb — suitable for LinkPrep/SV work

# Regions of interest — add or remove as needed
REGIONS = {
    "Sox11_Chr12":   "chr12:26,000,000-28,000,000",
    "Mir9-2_Chr13":  "chr13:83,500,000-84,500,000",
}

# Insulation score windows for TAD calling (bp)
INSULATION_WINDOWS = [50_000, 100_000, 200_000]

# SV detection: pairs with genomic distance above this are "long-range"
SV_DISTANCE_THRESHOLD = 2_000_000   # 2 Mb

# CNV: rolling-window size (number of bins) for smoothing
CNV_SMOOTHING_BINS = 5

# Phasing: minimum number of heterozygous SNPs per haplotype block
PHASING_MIN_SNPS = 3

# ---------------------------------------------------------------------------
# OUTPUT
# ---------------------------------------------------------------------------

FIGURES_DIR = os.path.join(os.path.dirname(__file__), "figures")
os.makedirs(FIGURES_DIR, exist_ok=True)

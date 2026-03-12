#!/usr/bin/env bash
# =============================================================================
# Central configuration for the HiChIP pipeline on Rutgers OARC Amarel.
# All SLURM scripts source this file.  Edit ONLY this file when adapting
# to new samples, directories, or resource allocations.
# =============================================================================

# ---------------------------------------------------------------------------
# YOUR IDENTITY — update before first run
# ---------------------------------------------------------------------------
NETID="your_netid"                              # Rutgers NetID
LAB_ACCOUNT="your_lab_account"                  # Slurm account (e.g. p_debs_1)
EMAIL="${NETID}@rutgers.edu"                    # Job notifications

# ---------------------------------------------------------------------------
# FILESYSTEM LAYOUT
# ---------------------------------------------------------------------------
# Amarel recommended layout:
#   /scratch/$NETID/  — large temporary storage (not backed up)
#   /projects/$LAB/   — project-level storage (backed up, quota applies)

SCRATCH="/scratch/${NETID}/hichip"              # working scratch space
PROJECT="/projects/${LAB_ACCOUNT}/hichip"       # long-term project storage
REPO_DIR="${PROJECT}/Dovetail-Analysis"         # cloned submodule scripts
REF_DIR="${PROJECT}/reference"                  # hg38 reference genome files
DATA_DIR="${SCRATCH}/fastq"                     # raw FASTQ input
BAM_DIR="${SCRATCH}/bam"                        # aligned BAMs
PAIRS_DIR="${SCRATCH}/pairs"                    # pairtools .pairs files
MATRIX_DIR="${SCRATCH}/matrices"                # cooler contact matrices
QC_DIR="${PROJECT}/qc"                          # QC reports
PEAKS_DIR="${PROJECT}/peaks"                    # MACS2 peak calls
LOOPS_DIR="${PROJECT}/loops"                    # FitHiChIP loop calls
LOGS_DIR="${PROJECT}/logs"                      # SLURM stdout/stderr

# ---------------------------------------------------------------------------
# REFERENCE GENOME (hg38)
# ---------------------------------------------------------------------------
REF_FASTA="${REF_DIR}/hg38.fa"
REF_GENOME="${REF_DIR}/hg38.genome"             # chrom sizes (2-col TSV)
BWA_INDEX="${REF_DIR}/hg38.fa"                  # bwa index prefix = fasta

# ---------------------------------------------------------------------------
# SAMPLES
# The pipeline is written as a SLURM array job.
# Add one entry per sample in the parallel arrays below.
# ---------------------------------------------------------------------------
SAMPLE_NAMES=(
    "CTCF_2M"          # small test dataset (~2M read pairs)
    "CTCF_DS"          # deep CTCF
    "H3K27Ac"          # active enhancer mark
    "H3K4me3"          # active promoter mark
)

# Paired-end FASTQs — must match SAMPLE_NAMES by index
R1_FASTQS=(
    "${DATA_DIR}/HiChiP_CTCF_2M_R1.fastq.gz"
    "${DATA_DIR}/CTCF-DS_R1.fastq.gz"
    "${DATA_DIR}/H3K27Ac_R1.fastq.gz"
    "${DATA_DIR}/H3K4me3_R1.fastq.gz"
)
R2_FASTQS=(
    "${DATA_DIR}/HiChiP_CTCF_2M_R2.fastq.gz"
    "${DATA_DIR}/CTCF-DS_R2.fastq.gz"
    "${DATA_DIR}/H3K27Ac_R2.fastq.gz"
    "${DATA_DIR}/H3K4me3_R2.fastq.gz"
)

# S3 download URLs (Dovetail public datasets — hg38, GM12878)
S3_BASE="https://s3.amazonaws.com/dovetail.pub/HiChIP/fastqs"
S3_URLS=(
    "${S3_BASE}/HiChiP_CTCF_2M_R1.fastq.gz"
    "${S3_BASE}/HiChiP_CTCF_2M_R2.fastq.gz"
    "${S3_BASE}/CTCF-DS_R1.fastq.gz"
    "${S3_BASE}/CTCF-DS_R2.fastq.gz"
    "${S3_BASE}/H3K27Ac_R1.fastq.gz"
    "${S3_BASE}/H3K27Ac_R2.fastq.gz"
    "${S3_BASE}/H3K4me3_R1.fastq.gz"
    "${S3_BASE}/H3K4me3_R2.fastq.gz"
)

# ENCODE ChIP-seq peak files for each sample (hg38, GM12878)
# These are used for FitHiChIP bias correction and ChIP enrichment QC.
# Set to "" for samples where you'll use MACS2 peaks from the BAM itself.
PEAK_FILES=(
    "${PEAKS_DIR}/CTCF_GM12878_ENCFF026ANI.bed"   # CTCF peaks
    "${PEAKS_DIR}/CTCF_GM12878_ENCFF026ANI.bed"   # CTCF peaks (same target)
    "${PEAKS_DIR}/H3K27Ac_GM12878_ENCFF254NLU.bed"
    "${PEAKS_DIR}/H3K4me3_GM12878_ENCFF659OHD.bed"
)

# FitHiChIP interaction type per sample:
#   1 = peak-to-peak  (CTCF, CTCF-binding-factor loops)
#   3 = peak-to-all   (H3K27ac, H3K4me3 — enhancer/promoter interactions)
FITHICHIP_INTTYPES=(1 1 3 3)

# ---------------------------------------------------------------------------
# SLURM DEFAULTS  (override per-script as needed)
# ---------------------------------------------------------------------------
PARTITION="main"
THREADS=16
MEM_ALIGN="64G"
MEM_PAIRS="32G"
MEM_MATRIX="64G"
MEM_LOOP="128G"
TIME_ALIGN="12:00:00"
TIME_PAIRS="8:00:00"
TIME_MATRIX="6:00:00"
TIME_LOOP="24:00:00"

# ---------------------------------------------------------------------------
# TOOL PATHS  (loaded via module or conda)
# ---------------------------------------------------------------------------
CONDA_ENV="hichip"                              # conda env name
SINGULARITY_IMG="${PROJECT}/containers/fithichip.sif"   # FitHiChIP container

# FitHiChIP analysis parameters
BINSIZE=2500
LOW_DIST=20000
UPP_DIST=2000000
QVALUE=0.01

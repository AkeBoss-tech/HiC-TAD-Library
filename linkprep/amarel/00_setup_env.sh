#!/usr/bin/env bash
# =============================================================================
# 00_setup_env.sh — One-time environment setup on Amarel
#
# Run this ONCE from a login node (not as a SLURM job — it's fast):
#   bash linkprep/amarel/00_setup_env.sh
#
# What it does:
#   1. Creates the conda environment "hichip" with all bioinformatics tools
#   2. Creates project/scratch directory structure
#   3. Pulls the FitHiChIP Singularity image
#   4. Downloads ENCODE peak files for GM12878
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

echo "=== HiChIP Environment Setup on Amarel ==="
echo "NetID    : ${NETID}"
echo "Scratch  : ${SCRATCH}"
echo "Project  : ${PROJECT}"
echo

# ---------------------------------------------------------------------------
# 1. Directory structure
# ---------------------------------------------------------------------------
echo "[1/5] Creating directory structure …"
mkdir -p "${DATA_DIR}" "${BAM_DIR}" "${PAIRS_DIR}" "${MATRIX_DIR}" \
         "${QC_DIR}" "${PEAKS_DIR}" "${LOOPS_DIR}" "${LOGS_DIR}" \
         "${REF_DIR}" "${PROJECT}/containers"
echo "  Done."

# ---------------------------------------------------------------------------
# 2. Conda environment
# ---------------------------------------------------------------------------
echo "[2/5] Creating conda environment '${CONDA_ENV}' …"
module purge
module load Anaconda3

conda create -y -n "${CONDA_ENV}" -c conda-forge -c bioconda \
    python=3.10 \
    bwa=0.7.17 \
    samtools=1.17 \
    pairtools=1.0.2 \
    cooler=0.9.3 \
    macs2=2.2.9.1 \
    deeptools=3.5.4 \
    preseq=3.2.0 \
    bedtools=2.31.0 \
    tabulate \
    numpy scipy pandas matplotlib \
    pysam pyBigWig

echo "  Conda env '${CONDA_ENV}' ready."

# ---------------------------------------------------------------------------
# 3. FitHiChIP Singularity container
# ---------------------------------------------------------------------------
echo "[3/5] Pulling FitHiChIP Singularity image …"
module load singularity
singularity pull --name "${SINGULARITY_IMG}" \
    docker://ay434/fithichip:latest
echo "  Image at: ${SINGULARITY_IMG}"

# ---------------------------------------------------------------------------
# 4. hg38 reference genome (skip if already present)
# ---------------------------------------------------------------------------
echo "[4/5] Checking reference genome …"
if [[ ! -f "${REF_FASTA}" ]]; then
    echo "  Downloading hg38 (UCSC) …"
    wget -q -O "${REF_FASTA}.gz" \
        https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
    gunzip "${REF_FASTA}.gz"
    samtools faidx "${REF_FASTA}"
    cut -f1,2 "${REF_FASTA}.fai" > "${REF_GENOME}"
    bwa index "${REF_FASTA}"
    echo "  hg38 indexed."
else
    echo "  hg38 already present — skipping."
fi

# ---------------------------------------------------------------------------
# 5. ENCODE peak files for GM12878 (hg38)
# ---------------------------------------------------------------------------
echo "[5/5] Downloading ENCODE GM12878 peak files …"

declare -A PEAKS=(
    ["CTCF_GM12878_ENCFF026ANI.bed"]="https://www.encodeproject.org/files/ENCFF026ANI/@@download/ENCFF026ANI.bed.gz"
    ["H3K27Ac_GM12878_ENCFF254NLU.bed"]="https://www.encodeproject.org/files/ENCFF254NLU/@@download/ENCFF254NLU.bed.gz"
    ["H3K4me3_GM12878_ENCFF659OHD.bed"]="https://www.encodeproject.org/files/ENCFF659OHD/@@download/ENCFF659OHD.bed.gz"
)

for fname in "${!PEAKS[@]}"; do
    dest="${PEAKS_DIR}/${fname}"
    if [[ ! -f "${dest}" ]]; then
        url="${PEAKS[${fname}]}"
        echo "  Downloading ${fname} …"
        wget -q -O "${dest}.gz" "${url}"
        gunzip "${dest}.gz"
        # Keep only first 3 columns (chr start end) for BED compatibility
        cut -f1,2,3 "${dest}" > "${dest}.tmp" && mv "${dest}.tmp" "${dest}"
    else
        echo "  ${fname} already present — skipping."
    fi
done

echo
echo "=== Setup complete ==="
echo "Next step: sbatch linkprep/amarel/01_download_data.slurm"

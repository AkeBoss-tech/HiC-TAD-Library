#!/usr/bin/env bash
#SBATCH --job-name=hichip_matrix
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/04_matrix_%A_%a.out
#SBATCH --error=__LOGS_DIR__/04_matrix_%A_%a.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 04_contact_matrix.slurm — Build cooler contact matrices
#
# Generates:
#   .cool  at 1 kb (fine) and 5 kb (standard)
#   .mcool with balanced (ICE-normalised) resolutions for HiGlass / our pipeline
#
# Tool: cooler (HDF5-based, works natively with this project's analysis code)
#
# Submit:
#   sbatch --dependency=afterok:<pairs_jobid> linkprep/amarel/04_contact_matrix.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
DEDUP_PAIRS="${PAIRS_DIR}/${SAMPLE}.dedup.pairs.gz"
THREADS=${SLURM_CPUS_PER_TASK:-16}

COOL_5KB="${MATRIX_DIR}/${SAMPLE}.5kb.cool"
COOL_1KB="${MATRIX_DIR}/${SAMPLE}.1kb.cool"
MCOOL="${MATRIX_DIR}/${SAMPLE}.mcool"

echo "=== Contact matrix: ${SAMPLE} ==="
date

module purge
module load Anaconda3
conda activate "${CONDA_ENV}"

mkdir -p "${MATRIX_DIR}"

# ---------------------------------------------------------------------------
# Index pairs for random access
# ---------------------------------------------------------------------------
echo "[1/4] Indexing pairs (pairix) …"
if [[ ! -f "${DEDUP_PAIRS}.px2" ]]; then
    pairix "${DEDUP_PAIRS}"
fi

# ---------------------------------------------------------------------------
# 5 kb cool (primary resolution for TAD/compartment analysis)
# ---------------------------------------------------------------------------
echo "[2/4] Building 5 kb cool …"
cooler cload pairix \
    -p "${THREADS}" \
    "${REF_GENOME}:5000" \
    "${DEDUP_PAIRS}" \
    "${COOL_5KB}"

# ---------------------------------------------------------------------------
# 1 kb cool (high-resolution loop calling)
# ---------------------------------------------------------------------------
echo "[3/4] Building 1 kb cool …"
cooler cload pairix \
    -p "${THREADS}" \
    "${REF_GENOME}:1000" \
    "${DEDUP_PAIRS}" \
    "${COOL_1KB}"

# ---------------------------------------------------------------------------
# Multi-resolution .mcool with ICE balancing
# Resolutions: 1k, 2k, 5k, 10k, 25k, 50k, 100k, 250k, 500k, 1M
# ---------------------------------------------------------------------------
echo "[4/4] Zoomifying to .mcool (ICE balance) …"
cooler zoomify \
    --balance \
    -p "${THREADS}" \
    --resolutions 1000,2000,5000,10000,25000,50000,100000,250000,500000,1000000 \
    -o "${MCOOL}" \
    "${COOL_5KB}"

echo "=== Contact matrix complete: ${SAMPLE} ==="
echo "  5 kb cool : ${COOL_5KB}"
echo "  1 kb cool : ${COOL_1KB}"
echo "  mcool     : ${MCOOL}"
ls -lh "${COOL_5KB}" "${COOL_1KB}" "${MCOOL}"
date

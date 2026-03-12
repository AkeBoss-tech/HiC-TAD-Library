#!/usr/bin/env bash
#SBATCH --job-name=hichip_pairs
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/03_pairtools_%A_%a.out
#SBATCH --error=__LOGS_DIR__/03_pairtools_%A_%a.err
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 03_pairtools.slurm — Parse, sort, dedup, and split aligned BAM
#
# Pipeline:
#   BAM → pairtools parse → pairtools sort → pairtools dedup
#       → pairtools split → .pairs + .bam (deduplicated)
#
# Key parameters (Dovetail recommended):
#   --min-mapq 40        : discard low-quality alignments
#   --walks-policy 5unique : handle multi-ligation events
#   --max-inter-align-gap 30 : tolerate small alignment gaps
#   --mark-dups          : flag PCR duplicates (for stats)
#
# Outputs per sample:
#   <sample>.mapped.pairs.gz   — all valid pairs (compressed)
#   <sample>.dedup.pairs.gz    — deduplicated valid pairs
#   <sample>.mapped.PT.bam     — deduplicated BAM (for ChIP QC)
#   <sample>.stats.txt         — pairtools statistics (input to get_qc.py)
#
# Submit:
#   sbatch --dependency=afterok:<align_jobid> linkprep/amarel/03_pairtools.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
IN_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
THREADS=${SLURM_CPUS_PER_TASK:-16}
TMP_DIR="${SCRATCH}/tmp/${SAMPLE}_pairs"

OUT_PAIRS="${PAIRS_DIR}/${SAMPLE}.mapped.pairs.gz"
DEDUP_PAIRS="${PAIRS_DIR}/${SAMPLE}.dedup.pairs.gz"
OUT_BAM="${BAM_DIR}/${SAMPLE}.mapped.PT.bam"
STATS="${QC_DIR}/${SAMPLE}.stats.txt"

echo "=== pairtools: ${SAMPLE} ==="
date

module purge
module load Anaconda3
conda activate "${CONDA_ENV}"

mkdir -p "${PAIRS_DIR}" "${TMP_DIR}" "${QC_DIR}"

# ---------------------------------------------------------------------------
# Step 1: Parse ligation junctions from BAM
# ---------------------------------------------------------------------------
echo "[1/4] pairtools parse …"
pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --nproc-in "${THREADS}" \
    --nproc-out "${THREADS}" \
    --chroms-path "${REF_GENOME}" \
    "${IN_BAM}" \
| pairtools sort \
    --nproc "${THREADS}" \
    --tmpdir "${TMP_DIR}" \
    --output "${OUT_PAIRS}"

echo "  Parsed pairs: ${OUT_PAIRS}"

# ---------------------------------------------------------------------------
# Step 2: Deduplicate + capture statistics
# ---------------------------------------------------------------------------
echo "[2/4] pairtools dedup …"
pairtools dedup \
    --nproc-in "${THREADS}" \
    --nproc-out "${THREADS}" \
    --mark-dups \
    --output-stats "${STATS}" \
    --output "${DEDUP_PAIRS}" \
    "${OUT_PAIRS}"

echo "  Dedup pairs  : ${DEDUP_PAIRS}"
echo "  Stats file   : ${STATS}"

# ---------------------------------------------------------------------------
# Step 3: Split into .pairs + .bam
# ---------------------------------------------------------------------------
echo "[3/4] pairtools split …"
pairtools split \
    --nproc-in "${THREADS}" \
    --nproc-out "${THREADS}" \
    --output-pairs "${DEDUP_PAIRS/.gz/.split.gz}" \
    --output-sam - \
    "${DEDUP_PAIRS}" \
| samtools sort -@ "${THREADS}" -T "${TMP_DIR}/bamsort" -o "${OUT_BAM}"
samtools index "${OUT_BAM}"

echo "  Deduplicated BAM: ${OUT_BAM}"

# ---------------------------------------------------------------------------
# Step 4: Print QC summary using Dovetail's get_qc.py
# ---------------------------------------------------------------------------
echo "[4/4] QC summary …"
python3 "${REPO_DIR}/get_qc.py" -p "${STATS}" \
    | tee "${QC_DIR}/${SAMPLE}.qc_report.txt"

echo "=== pairtools complete: ${SAMPLE} ==="
date
rm -rf "${TMP_DIR}"

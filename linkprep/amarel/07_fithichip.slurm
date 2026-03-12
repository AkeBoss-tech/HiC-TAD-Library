#!/usr/bin/env bash
#SBATCH --job-name=hichip_loops
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/07_fithichip_%A_%a.out
#SBATCH --error=__LOGS_DIR__/07_fithichip_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=128G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 07_fithichip.slurm — HiChIP loop calling with FitHiChIP
#
# FitHiChIP workflow (Dovetail docs):
#   1. Convert .pairs → HiC-Pro format (pairs column reorder + gzip)
#   2. Fill config template with sample-specific paths
#   3. Run FitHiChIP via Singularity container
#
# Key outputs (in ${LOOPS_DIR}/<sample>/):
#   *_interactions_FitHiC_Q0.01.bed              — all significant interactions
#   *_interactions_FitHiC_Q0.01_MergeNearContacts.bed  — merged (use this)
#   FitHiChIP_summary.html                       — diagnostic report
#
# QC threshold: >500 significant merged loops at 2.5 kb resolution
#               indicates a good-quality library
#
# Submit:
#   sbatch --dependency=afterok:<pairs_jobid>:<peaks_jobid> \
#          linkprep/amarel/07_fithichip.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
INTTYPE="${FITHICHIP_INTTYPES[$IDX]}"
DEDUP_PAIRS="${PAIRS_DIR}/${SAMPLE}.dedup.pairs.gz"
THREADS=${SLURM_CPUS_PER_TASK:-16}

# Prefer ENCODE peaks; fall back to MACS2 peaks
ENCODE_PEAK="${PEAK_FILES[$IDX]}"
MACS2_PEAK="${PEAKS_DIR}/${SAMPLE}.macs2_peaks.bed"
if [[ -f "${ENCODE_PEAK}" ]]; then
    PEAK_FILE="${ENCODE_PEAK}"
    echo "  Using ENCODE peaks: ${ENCODE_PEAK}"
elif [[ -f "${MACS2_PEAK}" ]]; then
    PEAK_FILE="${MACS2_PEAK}"
    echo "  Using MACS2 peaks : ${MACS2_PEAK}"
else
    echo "ERROR: No peak file found for ${SAMPLE}. Run 06_macs2_peaks.slurm first."
    exit 1
fi

HICPRO_PAIRS="${PAIRS_DIR}/${SAMPLE}.hicpro.valid.pairs.gz"
LOOP_OUT="${LOOPS_DIR}/${SAMPLE}"
CONFIG_FILE="${LOOP_OUT}/fithichip_config.txt"

echo "=== FitHiChIP loop calling: ${SAMPLE} ==="
echo "  Interaction type: ${INTTYPE}  (1=P2P, 3=P2A)"
echo "  Peak file: ${PEAK_FILE}"
echo "  Resolution: ${BINSIZE} bp"
date

module purge
module load Anaconda3 singularity
conda activate "${CONDA_ENV}"

mkdir -p "${LOOP_OUT}"

# ---------------------------------------------------------------------------
# Step 1: Convert pairs → HiC-Pro format
# Dovetail: "grep -v '#' *.pairs | awk '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}'"
# Column mapping from 4DN pairs to HiC-Pro valid pairs:
#   $1=readID $2=chr1 $3=pos1 $4=chr2 $5=pos2 $6=strand1 $7=strand2
# HiC-Pro format: readID chr1 pos1 strand1 chr2 pos2 strand2
# ---------------------------------------------------------------------------
echo "[1/3] Converting pairs to HiC-Pro format …"
zcat "${DEDUP_PAIRS}" \
    | grep -v '^#' \
    | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5"\t"$7}' \
    | gzip -c \
    > "${HICPRO_PAIRS}"
echo "  HiC-Pro pairs: ${HICPRO_PAIRS}"

# ---------------------------------------------------------------------------
# Step 2: Fill FitHiChIP config template
# ---------------------------------------------------------------------------
echo "[2/3] Building FitHiChIP config …"
sed \
    -e "s|__HICPRO_PAIRS__|${HICPRO_PAIRS}|g" \
    -e "s|__PEAK_FILE__|${PEAK_FILE}|g" \
    -e "s|__OUT_DIR__|${LOOP_OUT}|g" \
    -e "s|__GENOME_FILE__|${REF_GENOME}|g" \
    -e "s|__PREFIX__|${SAMPLE}.${BINSIZE}bp|g" \
    -e "s|__INTTYPE__|${INTTYPE}|g" \
    -e "s|__BINSIZE__|${BINSIZE}|g" \
    -e "s|__LOW_DIST__|${LOW_DIST}|g" \
    -e "s|__UPP_DIST__|${UPP_DIST}|g" \
    -e "s|__QVALUE__|${QVALUE}|g" \
    "$(dirname "$0")/fithichip_config.template.txt" \
    > "${CONFIG_FILE}"
echo "  Config: ${CONFIG_FILE}"

# ---------------------------------------------------------------------------
# Step 3: Run FitHiChIP via Singularity
# The container mounts:
#   /scratch, /projects → pass-through via $SINGULARITYENV
# FitHiChIP_Docker.sh is the entrypoint inside the container
# ---------------------------------------------------------------------------
echo "[3/3] Running FitHiChIP …"
singularity exec \
    --bind "${SCRATCH}:${SCRATCH}" \
    --bind "${PROJECT}:${PROJECT}" \
    "${SINGULARITY_IMG}" \
    bash -c "FitHiChIP_Docker.sh -C ${CONFIG_FILE}"

# ---------------------------------------------------------------------------
# Summary
# ---------------------------------------------------------------------------
MERGED_LOOPS="${LOOP_OUT}/${SAMPLE}.${BINSIZE}bp_interactions_FitHiC_Q${QVALUE}_MergeNearContacts.bed"
if [[ -f "${MERGED_LOOPS}" ]]; then
    N_LOOPS=$(wc -l < "${MERGED_LOOPS}")
    echo "=== FitHiChIP complete: ${SAMPLE} ==="
    echo "  Merged loops: ${N_LOOPS}  (file: ${MERGED_LOOPS})"
    if (( N_LOOPS > 500 )); then
        echo "  [PASS] Loop QC: >500 significant merged loops"
    else
        echo "  [WARN] Loop QC: only ${N_LOOPS} loops — consider deeper sequencing"
    fi
else
    echo "  WARNING: Expected output not found: ${MERGED_LOOPS}"
    echo "  Check FitHiChIP log in: ${LOOP_OUT}"
fi
date

#!/usr/bin/env bash
#SBATCH --job-name=hichip_diff
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/08_comparative_%j.out
#SBATCH --error=__LOGS_DIR__/08_comparative_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 08_comparative.slurm — Differential loop analysis between two conditions
#
# Identifies loops that are gained, lost, or unchanged between conditions.
# Uses FitHiChIP's built-in differential analysis module (EdgeR-based).
#
# Dovetail docs:
#   "Detect differentially interacting regions across conditions (treated/
#    untreated, WT/KO, normal/tumor)"
#   Output: Loops_EdgeR_Default_SIG.bed
#   Key columns: logFC (fold-change), FDR (adjusted p-value)
#
# USAGE:
#   Edit CONDITION_A and CONDITION_B below (must be sample names from config.sh),
#   then submit:
#   sbatch --dependency=afterok:<fithichip_jobid> linkprep/amarel/08_comparative.slurm
#
# NOTE: This script compares CTCF_2M vs CTCF_DS as a demonstration.
#       For your experiment, change A/B to e.g. "H3K27Ac_WT" vs "H3K27Ac_KO".
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

CONDITION_A="CTCF_2M"    # control / reference condition
CONDITION_B="CTCF_DS"    # treatment / comparison condition

LOOPS_A="${LOOPS_DIR}/${CONDITION_A}/${CONDITION_A}.${BINSIZE}bp_interactions_FitHiC.bed"
LOOPS_B="${LOOPS_DIR}/${CONDITION_B}/${CONDITION_B}.${BINSIZE}bp_interactions_FitHiC.bed"
DIFF_OUT="${LOOPS_DIR}/differential_${CONDITION_A}_vs_${CONDITION_B}"

echo "=== Differential loop analysis ==="
echo "  Condition A (reference) : ${CONDITION_A}"
echo "  Condition B (comparison): ${CONDITION_B}"
date

module purge
module load Anaconda3 singularity
conda activate "${CONDA_ENV}"

mkdir -p "${DIFF_OUT}"

# Validate inputs
for f in "${LOOPS_A}" "${LOOPS_B}"; do
    if [[ ! -f "$f" ]]; then
        echo "ERROR: Loop file not found: $f"
        echo "Run 07_fithichip.slurm first."
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# FitHiChIP differential analysis
# The diff module compares two loop files using EdgeR (count-based stats).
# For single-replicate comparisons, it uses the square-root dispersion method.
# ---------------------------------------------------------------------------
echo "[1/2] Running FitHiChIP differential analysis …"
singularity exec \
    --bind "${SCRATCH}:${SCRATCH}" \
    --bind "${PROJECT}:${PROJECT}" \
    "${SINGULARITY_IMG}" \
    bash -c "
        Rscript /FitHiChIP/src/Diff_Loops.r \
            --FileList ${LOOPS_A},${LOOPS_B} \
            --Labels ${CONDITION_A},${CONDITION_B} \
            --OutDir ${DIFF_OUT} \
            --FDRThr 0.05 \
            --logFCThr 1.0
    "

# ---------------------------------------------------------------------------
# Parse and summarise results
# ---------------------------------------------------------------------------
echo "[2/2] Summarising results …"
SIG_LOOPS="${DIFF_OUT}/Loops_EdgeR_Default_SIG.bed"
if [[ -f "${SIG_LOOPS}" ]]; then
    TOTAL=$(wc -l < "${SIG_LOOPS}")
    GAINED=$(awk '$NF > 0' "${SIG_LOOPS}" | wc -l)   # positive logFC
    LOST=$(awk '$NF < 0' "${SIG_LOOPS}" | wc -l)      # negative logFC
    echo
    echo "  Differential loops: ${TOTAL} total"
    echo "    Gained in ${CONDITION_B}: ${GAINED}"
    echo "    Lost   in ${CONDITION_B}: ${LOST}"
    echo "  Output: ${SIG_LOOPS}"
else
    echo "  WARNING: Significant loops file not found — check FitHiChIP diff log."
fi

echo "=== Differential analysis complete ==="
date

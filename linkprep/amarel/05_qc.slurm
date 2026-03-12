#!/usr/bin/env bash
#SBATCH --job-name=hichip_qc
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/05_qc_%A_%a.out
#SBATCH --error=__LOGS_DIR__/05_qc_%A_%a.err
#SBATCH --time=02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 05_qc.slurm — HiChIP quality control assessment
#
# Runs three QC layers:
#
#   Layer 1 — Proximity-ligation QC (from pairtools stats)
#     Tool : Dovetail get_qc.py
#     Pass : No-dup cis ≥1kb > 20% of total no-dup pairs
#
#   Layer 2 — ChIP enrichment QC (reads in peak windows)
#     Tool : Dovetail enrichment_stats.sh
#     Pass : >2% reads in 1000 bp window around peaks
#
#   Layer 3 — Library complexity projection
#     Tool : preseq lc_extrap
#     Tells you how many unique pairs you'd get at 400M+ read depth
#
# Submit:
#   sbatch --dependency=afterok:<pairs_jobid> linkprep/amarel/05_qc.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
STATS="${QC_DIR}/${SAMPLE}.stats.txt"
COORD_BAM="${BAM_DIR}/${SAMPLE}.coord.bam"
PEAK_FILE="${PEAK_FILES[$IDX]}"
THREADS=${SLURM_CPUS_PER_TASK:-16}
QC_OUT="${QC_DIR}/${SAMPLE}"

echo "=== QC: ${SAMPLE} ==="
date

module purge
module load Anaconda3
conda activate "${CONDA_ENV}"

mkdir -p "${QC_DIR}"

# ---------------------------------------------------------------------------
# Layer 1: Proximity-ligation stats (get_qc.py from Dovetail repo)
# ---------------------------------------------------------------------------
echo "[1/4] Proximity-ligation QC …"
python3 "${REPO_DIR}/get_qc.py" -p "${STATS}" \
    | tee "${QC_OUT}.pairligation_qc.txt"

# Parse pass/fail thresholds
CIS_1KB=$(grep "cis_1kb+" "${STATS}" | awk '{print $2}')
TOTAL_NODUPS=$(grep "^total_nodups" "${STATS}" | awk '{print $2}')
if [[ -n "${CIS_1KB}" && -n "${TOTAL_NODUPS}" && "${TOTAL_NODUPS}" -gt 0 ]]; then
    PCT=$(echo "scale=2; ${CIS_1KB}*100/${TOTAL_NODUPS}" | bc)
    echo "  No-dup cis >=1kb: ${PCT}%  (threshold: >20%)"
    if (( $(echo "${PCT} > 20" | bc -l) )); then
        echo "  [PASS] Proximity-ligation QC"
    else
        echo "  [FAIL] Proximity-ligation QC — consider re-sequencing"
    fi
fi

# ---------------------------------------------------------------------------
# Layer 2: ChIP enrichment (enrichment_stats.sh from Dovetail repo)
# ---------------------------------------------------------------------------
echo "[2/4] ChIP enrichment QC …"
ENRICH_SCRIPT="${REPO_DIR}/docs/source/hichip/enrichment_stats.sh"
if [[ -f "${ENRICH_SCRIPT}" && -f "${PEAK_FILE}" ]]; then
    chmod +x "${ENRICH_SCRIPT}"
    "${ENRICH_SCRIPT}" \
        -g "${REF_GENOME}" \
        -b "${COORD_BAM}" \
        -p "${PEAK_FILE}" \
        -t "${THREADS}" \
        -x "${QC_OUT}" \
        | tee "${QC_OUT}.enrichment_qc.txt"
else
    echo "  enrichment_stats.sh or peak file not found — skipping."
    echo "  (peak file: ${PEAK_FILE})"
fi

# ---------------------------------------------------------------------------
# Layer 3: Library complexity projection (preseq)
# ---------------------------------------------------------------------------
echo "[3/4] Library complexity projection (preseq) …"
# preseq needs sorted, duplicate-marked BAM
# Use the coordinate-sorted BAM from step 02
preseq lc_extrap \
    -bam \
    -pe \
    -output "${QC_OUT}.preseq.txt" \
    "${COORD_BAM}" \
    2> "${QC_OUT}.preseq.log" || echo "  preseq warning (see log)"

# ---------------------------------------------------------------------------
# Layer 4: ChIP enrichment plot (Dovetail script)
# ---------------------------------------------------------------------------
echo "[4/4] ChIP enrichment plot …"
PLOT_SCRIPT="${REPO_DIR}/docs/source/hichip/plot_chip_enrichment_bed.py"
if [[ -f "${PLOT_SCRIPT}" && -f "${PEAK_FILE}" ]]; then
    python3 "${PLOT_SCRIPT}" \
        -bam "${COORD_BAM}" \
        -peaks "${PEAK_FILE}" \
        -output "${QC_OUT}.enrichment.png"
    echo "  Plot: ${QC_OUT}.enrichment.png"
fi

echo "=== QC complete: ${SAMPLE} ==="
echo "  Reports in: ${QC_DIR}"
ls -lh "${QC_DIR}/${SAMPLE}"* 2>/dev/null || true
date

#!/usr/bin/env bash
# =============================================================================
# submit_pipeline.sh — Submit the full HiChIP pipeline to Amarel SLURM
#
# Uses SLURM job dependencies so each step only starts after the previous
# step completes successfully (afterok dependency chain).
#
# Usage:
#   # Dry run — print job commands without submitting
#   bash linkprep/amarel/submit_pipeline.sh --dry-run
#
#   # Full run (starts with data download)
#   bash linkprep/amarel/submit_pipeline.sh
#
#   # Skip download (data already present)
#   bash linkprep/amarel/submit_pipeline.sh --skip-download
#
#   # Run specific steps only (comma-separated)
#   bash linkprep/amarel/submit_pipeline.sh --steps 02,03,04
#
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

AMAREL_DIR="$(dirname "$0")"
DRY_RUN=false
SKIP_DOWNLOAD=false
STEPS="01,02,03,04,05,06,07"   # 08 (comparative) is run separately

# Parse arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --dry-run)       DRY_RUN=true ;;
        --skip-download) SKIP_DOWNLOAD=true ;;
        --steps)         STEPS="$2"; shift ;;
        *) echo "Unknown argument: $1"; exit 1 ;;
    esac
    shift
done

# Determine array bounds from number of samples
N_SAMPLES=${#SAMPLE_NAMES[@]}
ARRAY_RANGE="0-$((N_SAMPLES - 1))"

# Replace config.sh placeholders in SLURM headers
prepare_script() {
    local src="$1"
    local tmp="${LOGS_DIR}/$(basename "$src").prepared"
    mkdir -p "${LOGS_DIR}"
    sed \
        -e "s|__LAB_ACCOUNT__|${LAB_ACCOUNT}|g" \
        -e "s|__LOGS_DIR__|${LOGS_DIR}|g" \
        -e "s|__EMAIL__|${EMAIL}|g" \
        "${src}" > "${tmp}"
    echo "${tmp}"
}

submit() {
    local script="$1"
    local dep="$2"
    local prepared
    prepared=$(prepare_script "${script}")

    local dep_flag=""
    if [[ -n "${dep}" ]]; then
        dep_flag="--dependency=afterok:${dep}"
    fi

    if $DRY_RUN; then
        echo "[DRY RUN] sbatch ${dep_flag} ${prepared}"
        echo "fake_jobid_$(basename "$script")"
    else
        sbatch ${dep_flag} "${prepared}" | awk '{print $NF}'
    fi
}

echo "=== HiChIP Pipeline Submission ==="
echo "  Samples       : ${SAMPLE_NAMES[*]}"
echo "  Array range   : ${ARRAY_RANGE}"
echo "  Steps         : ${STEPS}"
echo "  Dry run       : ${DRY_RUN}"
echo

JOB_PREV=""

# Step 01 — Download
if echo "${STEPS}" | grep -q "01" && ! $SKIP_DOWNLOAD; then
    echo "[01] Submitting data download …"
    JOB_DL=$(submit "${AMAREL_DIR}/01_download_data.slurm" "")
    echo "  Job ID: ${JOB_DL}"
    JOB_PREV="${JOB_DL}"
fi

# Step 02 — Align (array)
if echo "${STEPS}" | grep -q "02"; then
    echo "[02] Submitting alignment (array ${ARRAY_RANGE}) …"
    JOB_ALIGN=$(submit "${AMAREL_DIR}/02_align.slurm" "${JOB_PREV}")
    echo "  Job ID: ${JOB_ALIGN}"
    JOB_PREV="${JOB_ALIGN}"
fi

# Step 03 — pairtools (array)
if echo "${STEPS}" | grep -q "03"; then
    echo "[03] Submitting pairtools (array ${ARRAY_RANGE}) …"
    JOB_PAIRS=$(submit "${AMAREL_DIR}/03_pairtools.slurm" "${JOB_PREV}")
    echo "  Job ID: ${JOB_PAIRS}"
    JOB_PREV="${JOB_PAIRS}"
fi

# Steps 04, 05, 06 can run in parallel after pairtools
PAIRS_JOB="${JOB_PREV}"

if echo "${STEPS}" | grep -q "04"; then
    echo "[04] Submitting contact matrix (array ${ARRAY_RANGE}) …"
    JOB_MATRIX=$(submit "${AMAREL_DIR}/04_contact_matrix.slurm" "${PAIRS_JOB}")
    echo "  Job ID: ${JOB_MATRIX}"
fi

if echo "${STEPS}" | grep -q "05"; then
    echo "[05] Submitting QC (array ${ARRAY_RANGE}) …"
    JOB_QC=$(submit "${AMAREL_DIR}/05_qc.slurm" "${PAIRS_JOB}")
    echo "  Job ID: ${JOB_QC}"
fi

if echo "${STEPS}" | grep -q "06"; then
    echo "[06] Submitting MACS2 peaks (array ${ARRAY_RANGE}) …"
    JOB_PEAKS=$(submit "${AMAREL_DIR}/06_macs2_peaks.slurm" "${PAIRS_JOB}")
    echo "  Job ID: ${JOB_PEAKS}"
    JOB_PREV="${JOB_PEAKS}"   # FitHiChIP needs peaks
fi

# Step 07 — FitHiChIP (depends on pairs + peaks)
if echo "${STEPS}" | grep -q "07"; then
    DEP_07="${PAIRS_JOB}"
    if [[ -n "${JOB_PEAKS:-}" ]]; then
        DEP_07="${PAIRS_JOB}:${JOB_PEAKS}"
    fi
    echo "[07] Submitting FitHiChIP loop calling (array ${ARRAY_RANGE}) …"
    JOB_LOOPS=$(submit "${AMAREL_DIR}/07_fithichip.slurm" "${DEP_07}")
    echo "  Job ID: ${JOB_LOOPS}"
fi

echo
echo "=== All jobs submitted ==="
echo "Monitor with:"
echo "  squeue -u ${NETID}"
echo "  sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS"
echo
echo "Comparative analysis (step 08) — run separately after step 07:"
echo "  sbatch --dependency=afterok:<fithichip_jobid> linkprep/amarel/08_comparative.slurm"

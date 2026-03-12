#!/usr/bin/env bash
#SBATCH --job-name=hichip_peaks
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/06_macs2_%A_%a.out
#SBATCH --error=__LOGS_DIR__/06_macs2_%A_%a.err
#SBATCH --time=04:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=0-3
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 06_macs2_peaks.slurm — 1D peak calling from HiChIP BAM using MACS2
#
# Context (from Dovetail docs):
#   ChIP-seq peaks are preferable over de novo MACS2 calls when available.
#   When ChIP-seq data is absent, extract primary alignments from the HiChIP
#   BAM, convert to BED, and run MACS2.
#
#   Caveats:
#   - May call secondary/artefact peaks without ChIP reference
#   - Run FitHiChIP with BOTH peak-to-peak AND peak-to-all to cross-check
#
# Outputs:
#   <sample>_peaks.narrowPeak   — MACS2 peak calls
#   <sample>_peaks.bed          — simplified 3-column BED (for FitHiChIP)
#
# Submit:
#   sbatch --dependency=afterok:<pairs_jobid> linkprep/amarel/06_macs2_peaks.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
COORD_BAM="${BAM_DIR}/${SAMPLE}.coord.bam"
THREADS=${SLURM_CPUS_PER_TASK:-8}
MACS2_OUT="${PEAKS_DIR}/${SAMPLE}_macs2"
PRIMARY_BED="${PEAKS_DIR}/${SAMPLE}.primary.aln.bed"

echo "=== MACS2 peak calling: ${SAMPLE} ==="
date

module purge
module load Anaconda3
conda activate "${CONDA_ENV}"

mkdir -p "${PEAKS_DIR}"

# ---------------------------------------------------------------------------
# Step 1: Extract primary alignments only
# Dovetail: "samtools -view -h -F 0x900 mapped.bam | bedtools bamtobed"
#   Flag 0x900 = SUPPLEMENTARY (0x800) + SECONDARY (0x100)
# ---------------------------------------------------------------------------
echo "[1/3] Extracting primary alignments …"
samtools view -h -F 0x900 "${COORD_BAM}" \
    | bedtools bamtobed -i stdin \
    > "${PRIMARY_BED}"
echo "  Primary alignments: $(wc -l < "$PRIMARY_BED") reads"

# ---------------------------------------------------------------------------
# Step 2: MACS2 peak calling
# ---------------------------------------------------------------------------
echo "[2/3] MACS2 callpeak …"
macs2 callpeak \
    -t "${PRIMARY_BED}" \
    -n "${SAMPLE}" \
    --outdir "${MACS2_OUT}" \
    -g hs \
    --nomodel \
    --extsize 147 \
    2> "${MACS2_OUT}/macs2.log"

# ---------------------------------------------------------------------------
# Step 3: Convert narrowPeak → simple 3-col BED for FitHiChIP
# ---------------------------------------------------------------------------
echo "[3/3] Converting to BED …"
NARROW="${MACS2_OUT}/${SAMPLE}_peaks.narrowPeak"
SIMPLE_BED="${PEAKS_DIR}/${SAMPLE}.macs2_peaks.bed"
if [[ -f "${NARROW}" ]]; then
    cut -f1,2,3 "${NARROW}" > "${SIMPLE_BED}"
    echo "  Peak BED: ${SIMPLE_BED}  ($(wc -l < "$SIMPLE_BED") peaks)"
else
    echo "  WARNING: narrowPeak not found — MACS2 may have failed (see log)."
fi

echo "=== MACS2 complete: ${SAMPLE} ==="
echo "  Output dir: ${MACS2_OUT}"
date

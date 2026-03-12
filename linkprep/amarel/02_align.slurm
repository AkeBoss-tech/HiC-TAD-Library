#!/usr/bin/env bash
#SBATCH --job-name=hichip_align
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/02_align_%A_%a.out
#SBATCH --error=__LOGS_DIR__/02_align_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --array=0-3             # one task per sample (CTCF_2M, CTCF_DS, H3K27Ac, H3K4me3)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 02_align.slurm — BWA-MEM alignment for HiChIP data
#
# Uses HiChIP-specific flags:
#   -5  report supplementary alignments for split reads
#   -S  skip mate rescue
#   -P  skip pairing; treat mates independently
#   -T0 report all alignments (required for proximity ligation reads)
#
# Output: sorted, indexed BAM per sample
#
# Submit after download:
#   sbatch --dependency=afterok:<download_jobid> linkprep/amarel/02_align.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

# SLURM array index → sample
IDX=${SLURM_ARRAY_TASK_ID}
SAMPLE="${SAMPLE_NAMES[$IDX]}"
R1="${R1_FASTQS[$IDX]}"
R2="${R2_FASTQS[$IDX]}"
OUT_BAM="${BAM_DIR}/${SAMPLE}.bam"
SORT_BAM="${BAM_DIR}/${SAMPLE}.sorted.bam"
THREADS=${SLURM_CPUS_PER_TASK:-16}

echo "=== Alignment: ${SAMPLE} ==="
echo "  R1  : ${R1}"
echo "  R2  : ${R2}"
echo "  Out : ${SORT_BAM}"
date

module purge
module load Anaconda3
conda activate "${CONDA_ENV}"

mkdir -p "${BAM_DIR}"

# Step 1: Align with BWA-MEM (HiChIP flags: -5SP -T0)
echo "[1/3] BWA-MEM alignment …"
bwa mem \
    -5SP \
    -T0 \
    -t "${THREADS}" \
    "${BWA_INDEX}" \
    "${R1}" \
    "${R2}" \
| samtools view -bS -F 4 - \
> "${OUT_BAM}"

# Step 2: Sort by name (required for pairtools parse)
echo "[2/3] Name-sorting BAM …"
samtools sort \
    -n \
    -@ "${THREADS}" \
    -T "${SCRATCH}/tmp/${SAMPLE}_sort" \
    -o "${SORT_BAM}" \
    "${OUT_BAM}"

rm -f "${OUT_BAM}"   # discard unsorted

# Step 3: Coordinate-sorted + indexed copy for ChIP QC (enrichment_stats.sh)
COORD_BAM="${BAM_DIR}/${SAMPLE}.coord.bam"
echo "[3/3] Coordinate-sorting for QC …"
samtools sort \
    -@ "${THREADS}" \
    -T "${SCRATCH}/tmp/${SAMPLE}_coord" \
    -o "${COORD_BAM}" \
    "${SORT_BAM}"
samtools index "${COORD_BAM}"

echo "=== Alignment complete: ${SAMPLE} ==="
echo "  Name-sorted : ${SORT_BAM}  ($(du -sh "$SORT_BAM" | cut -f1))"
echo "  Coord-sorted: ${COORD_BAM}  ($(du -sh "$COORD_BAM" | cut -f1))"
date

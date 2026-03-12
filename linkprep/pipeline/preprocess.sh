#!/usr/bin/env bash
# =============================================================================
# LinkPrep pre-processing pipeline: FASTQ → contact matrix (.cool)
#
# Usage:
#   bash pipeline/preprocess.sh \
#       <reference.fa> \
#       <read1.fastq.gz> \
#       <read2.fastq.gz> \
#       <sample_name> \
#       [threads=8]
#
# Requirements (conda install -c bioconda):
#   bwa samtools pairtools cooler
#
# When real data arrives, just call this script with your paths.
# =============================================================================

set -euo pipefail

REF=${1:?"Usage: $0 <ref.fa> <R1.fq.gz> <R2.fq.gz> <sample> [threads]"}
R1=${2:?"Usage: $0 <ref.fa> <R1.fq.gz> <R2.fq.gz> <sample> [threads]"}
R2=${3:?"Usage: $0 <ref.fa> <R1.fq.gz> <R2.fq.gz> <sample> [threads]"}
SAMPLE=${4:?"Usage: $0 <ref.fa> <R1.fq.gz> <R2.fq.gz> <sample> [threads]"}
THREADS=${5:-8}

OUTDIR="sample_data/output"
mkdir -p "$OUTDIR"

BIN_SIZE=10000   # 10 kb — adjust to 5000 for finer TAD resolution

echo "=== LinkPrep Pre-processing: $SAMPLE ==="
echo "Reference : $REF"
echo "Read 1    : $R1"
echo "Read 2    : $R2"
echo "Threads   : $THREADS"
echo "Bin size  : $BIN_SIZE bp"
echo

# ---------------------------------------------------------------------------
# Step 1 – Index reference (skip if already done)
# ---------------------------------------------------------------------------
if [[ ! -f "${REF}.bwt" ]]; then
    echo "[1/6] Indexing reference genome …"
    bwa index "$REF"
fi

# ---------------------------------------------------------------------------
# Step 2 – Align with BWA-MEM (LinkPrep-specific flags)
#   -5  : always use mate-rescue for chimeric reads
#   -S  : skip mate-rescue for supplementary alignments
#   -P  : skip pairing; output best alignment regardless of pairing
#   -T0 : report all hits (zero minimum alignment score)
# ---------------------------------------------------------------------------
echo "[2/6] Aligning reads (bwa mem -5SP -T0) …"
BAM="${OUTDIR}/${SAMPLE}.bam"
bwa mem -5SP -T0 -t "$THREADS" "$REF" "$R1" "$R2" \
    | samtools view -bS -F 4 - \
    | samtools sort -n -T "${OUTDIR}/${SAMPLE}_tmp" -o "$BAM"
echo "  → $BAM"

# ---------------------------------------------------------------------------
# Step 3 – Parse ligation junctions with pairtools
# ---------------------------------------------------------------------------
echo "[3/6] Parsing ligation junctions (pairtools parse) …"
CHROM_SIZES="${OUTDIR}/${SAMPLE}.chrom.sizes"
samtools view -H "$BAM" \
    | awk '/^@SQ/{split($2,a,":");split($3,b,":");print a[2]"\t"b[2]}' \
    > "$CHROM_SIZES"

PARSED="${OUTDIR}/${SAMPLE}.parsed.pairs.gz"
pairtools parse \
    --min-mapq 40 \
    --walks-policy 5unique \
    --max-inter-align-gap 30 \
    --nproc-in "$THREADS" \
    --nproc-out "$THREADS" \
    --chroms-path "$CHROM_SIZES" \
    "$BAM" \
    | pairtools sort \
        --nproc "$THREADS" \
        --tmpdir "${OUTDIR}" \
    | gzip > "$PARSED"
echo "  → $PARSED"

# ---------------------------------------------------------------------------
# Step 4 – Mark & remove PCR duplicates
# ---------------------------------------------------------------------------
echo "[4/6] Deduplicating (pairtools dedup) …"
DEDUP="${OUTDIR}/${SAMPLE}.dedup.pairs.gz"
STATS="${OUTDIR}/${SAMPLE}.dedup.stats.txt"
pairtools dedup \
    --nproc-in "$THREADS" \
    --nproc-out "$THREADS" \
    --mark-dups \
    --output-stats "$STATS" \
    --output "$DEDUP" \
    "$PARSED"
echo "  → $DEDUP"
echo "  Stats → $STATS"

# ---------------------------------------------------------------------------
# Step 5 – Keep only valid pairs (UU = both mates uniquely mapped)
# ---------------------------------------------------------------------------
echo "[5/6] Filtering valid pairs …"
VALID="${OUTDIR}/${SAMPLE}.valid.pairs.gz"
pairtools select '(pair_type == "UU") and (chrom1 == chrom2 or True)' \
    --output "$VALID" \
    "$DEDUP"

# Symlink as the canonical pairs file expected by config.py
CANONICAL="${OUTDIR}/linkprep_sample.pairs.gz"
ln -sf "$(basename "$VALID")" "$CANONICAL" 2>/dev/null || cp "$VALID" "$CANONICAL"
echo "  → $VALID"

# ---------------------------------------------------------------------------
# Step 6 – Generate contact matrix with Cooler
# ---------------------------------------------------------------------------
echo "[6/6] Building contact matrix (cooler cload pairs) …"
COOL="${OUTDIR}/linkprep_sample.cool"
cooler cload pairs \
    --chrom1 2 --pos1 3 --chrom2 4 --pos2 5 \
    "${CHROM_SIZES}:${BIN_SIZE}" \
    "$VALID" \
    "$COOL"

# Balance the matrix (ICE normalisation)
cooler balance --force "$COOL"

echo "  → $COOL"
echo
echo "=== Pre-processing complete ==="
echo "QC stats : $STATS"
echo "Cool file: $COOL"
echo
echo "Next steps:"
echo "  python pipeline/qc_metrics.py"
echo "  python run_analysis.py"

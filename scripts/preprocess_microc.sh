#!/usr/bin/env bash
# =============================================================================
# preprocess_microc.sh — Dovetail LinkPrep Micro-C preprocessing pipeline
#
# Converts raw paired FASTQ files into an ICE-normalised multi-resolution
# cooler file (.mcool) compatible with this library.
#
# Pipeline
# --------
#   FASTQ (R1, R2)
#     └─ bwa mem → SAM
#          └─ pairtools parse  → PAIRSAM     (ligation junction detection)
#               └─ pairtools sort            (coordinate sort)
#                    └─ pairtools dedup      (remove PCR duplicates)
#                         └─ pairtools split → .pairs + unsorted BAM
#                              ├─ samtools sort+index → mapped.PT.bam
#                              └─ bgzip + pairix → mapped.pairs.gz
#                                   └─ cooler cload + zoomify → output.mcool
#
# QC thresholds to check after running
# -------------------------------------
#   - No-dup cis pairs ≥ 1 kb  > 40 % of total no-dup pairs  (ligation quality)
#   - No-dup pairs              > 125 M for TAD/loop analysis
#   - Library complexity        > 125 M unique pairs at 400 M reads
#
# Usage
# -----
#   bash scripts/preprocess_microc.sh
#   bash scripts/preprocess_microc.sh --help
#
# Requirements (must be in PATH or conda env)
#   bwa, samtools, pairtools, bgzip, pairix, cooler
#
# Install in hic-analysis conda environment:
#   conda activate hic-analysis
#   conda install -c bioconda bwa samtools pairtools pairix cooler htslib
# =============================================================================

set -euo pipefail

# =============================================================================
# ── CONFIGURATION — edit these before running ────────────────────────────────
# =============================================================================

# Input reads (paired-end FASTQ, may be gzipped)
R1_FASTQ="${R1_FASTQ:-data/raw/sample_R1.fastq.gz}"
R2_FASTQ="${R2_FASTQ:-data/raw/sample_R2.fastq.gz}"

# Reference genome FASTA (must already be indexed with bwa index and samtools faidx)
GENOME_FA="${GENOME_FA:-data/raw/mm10.fa}"

# Chromosome sizes file  (generate with: samtools faidx mm10.fa; cut -f1,2 mm10.fa.fai > mm10.chromsizes)
CHROMSIZES="${CHROMSIZES:-data/raw/mm10.chromsizes}"

# Output directory
OUTPUT_DIR="${OUTPUT_DIR:-data/processed/microc_$(date +%Y%m%d)}"

# Sample name (used for output filenames)
SAMPLE="${SAMPLE:-sample}"

# Resolution for cooler (bp). Use 1000 for loops, 5000 for TADs.
# cooler zoomify will generate all coarser resolutions automatically.
BIN_SIZE="${BIN_SIZE:-1000}"

# Number of CPU threads for alignment and pairtools
CORES="${CORES:-$(nproc 2>/dev/null || sysctl -n hw.logicalcpu 2>/dev/null || echo 8)}"

# Minimum MAPQ to keep a read pair (40 = recommended for Micro-C)
MIN_MAPQ="${MIN_MAPQ:-40}"

# Temporary directory for pairtools sort (should have plenty of disk space)
TMPDIR="${TMPDIR:-/tmp}"

# =============================================================================
# ── HELPERS ──────────────────────────────────────────────────────────────────
# =============================================================================

log() { echo "[$(date '+%H:%M:%S')] $*"; }

check_tools() {
    local missing=()
    for tool in bwa samtools pairtools bgzip pairix cooler; do
        command -v "$tool" &>/dev/null || missing+=("$tool")
    done
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "ERROR: missing tools: ${missing[*]}"
        echo "  Install with: conda install -c bioconda ${missing[*]}"
        exit 1
    fi
}

print_help() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \?//'
    exit 0
}

check_inputs() {
    local missing=()
    [[ -f "$R1_FASTQ" ]] || missing+=("$R1_FASTQ")
    [[ -f "$R2_FASTQ" ]] || missing+=("$R2_FASTQ")
    [[ -f "$GENOME_FA" ]] || missing+=("$GENOME_FA")
    [[ -f "$CHROMSIZES" ]] || missing+=("$CHROMSIZES")
    if [[ ${#missing[@]} -gt 0 ]]; then
        echo "ERROR: missing input files:"
        printf '  %s\n' "${missing[@]}"
        exit 1
    fi
}

# =============================================================================
# ── MAIN ─────────────────────────────────────────────────────────────────────
# =============================================================================

[[ "${1:-}" == "--help" ]] && print_help

check_tools
check_inputs

mkdir -p "$OUTPUT_DIR"

PAIRSAM_RAW="$OUTPUT_DIR/${SAMPLE}.raw.pairsam.gz"
PAIRSAM_SORT="$OUTPUT_DIR/${SAMPLE}.sorted.pairsam.gz"
PAIRSAM_DEDUP="$OUTPUT_DIR/${SAMPLE}.dedup.pairsam.gz"
PAIRS_GZ="$OUTPUT_DIR/${SAMPLE}.mapped.pairs.gz"
BAM_UNSORTED="$OUTPUT_DIR/${SAMPLE}.unsorted.bam"
BAM_FINAL="$OUTPUT_DIR/${SAMPLE}.mapped.PT.bam"
DEDUP_STATS="$OUTPUT_DIR/${SAMPLE}.dedup.stats"
COOL_OUT="$OUTPUT_DIR/${SAMPLE}.${BIN_SIZE}.cool"
MCOOL_OUT="$OUTPUT_DIR/${SAMPLE}.mcool"

# ---------------------------------------------------------------------------
# Step 1: Alignment (bwa mem)
# Key flags:
#   -5  Use 5'-most alignment as primary (required for Hi-C/Micro-C)
#   -S  Skip mate rescue
#   -P  Skip mate pairing
#   -T0 Output all alignments regardless of score (pairtools will filter)
# ---------------------------------------------------------------------------
log "Step 1/7: Aligning with bwa mem (${CORES} threads)..."
bwa mem -5SP -T0 -t "$CORES" \
    "$GENOME_FA" "$R1_FASTQ" "$R2_FASTQ" \
    | pairtools parse \
        --min-mapq "$MIN_MAPQ" \
        --walks-policy 5unique \
        --max-inter-align-gap 30 \
        --nproc-in "$CORES" --nproc-out "$CORES" \
        --chroms-path "$CHROMSIZES" \
        -o "$PAIRSAM_RAW"

log "  → $PAIRSAM_RAW"

# ---------------------------------------------------------------------------
# Step 2: Sort PAIRSAM by genomic coordinates
# Note: sorting is I/O-heavy; ensure TMPDIR has sufficient disk space
# ---------------------------------------------------------------------------
log "Step 2/7: Sorting PAIRSAM..."
pairtools sort \
    --nproc "$CORES" \
    --tmpdir "$TMPDIR" \
    -o "$PAIRSAM_SORT" \
    "$PAIRSAM_RAW"

rm -f "$PAIRSAM_RAW"  # free disk space
log "  → $PAIRSAM_SORT"

# ---------------------------------------------------------------------------
# Step 3: Remove PCR duplicates
# Duplicates are identified by identical 5'-end coordinates on both sides.
# They are marked as pair-type "DD" and written to a separate .dups file.
# ---------------------------------------------------------------------------
log "Step 3/7: Removing PCR duplicates..."
pairtools dedup \
    --nproc-in "$CORES" --nproc-out "$CORES" \
    --mark-dups \
    --output-stats "$DEDUP_STATS" \
    -o "$PAIRSAM_DEDUP" \
    "$PAIRSAM_SORT"

rm -f "$PAIRSAM_SORT"
log "  → $PAIRSAM_DEDUP  (stats: $DEDUP_STATS)"

# QC check: print duplicate rate
if [[ -f "$DEDUP_STATS" ]]; then
    TOTAL=$(grep '^total ' "$DEDUP_STATS" | awk '{print $2}')
    DUPS=$(grep '^duplicate ' "$DEDUP_STATS" | awk '{print $2}')
    if [[ -n "$TOTAL" && "$TOTAL" -gt 0 ]]; then
        DUP_PCT=$(awk "BEGIN{printf \"%.1f\", $DUPS/$TOTAL*100}")
        log "  Duplicate rate: ${DUP_PCT}%  (${DUPS}/${TOTAL})"
    fi
fi

# ---------------------------------------------------------------------------
# Step 4: Split into .pairs and BAM
# ---------------------------------------------------------------------------
log "Step 4/7: Splitting to .pairs and BAM..."
pairtools split \
    --nproc-in "$CORES" --nproc-out "$CORES" \
    --output-pairs "$PAIRS_GZ" \
    --output-sam "$BAM_UNSORTED" \
    "$PAIRSAM_DEDUP"

rm -f "$PAIRSAM_DEDUP"
log "  → $PAIRS_GZ"
log "  → $BAM_UNSORTED"

# ---------------------------------------------------------------------------
# Step 5: Sort and index BAM
# Produces mapped.PT.bam for any downstream variant calling or visualization
# ---------------------------------------------------------------------------
log "Step 5/7: Sorting and indexing BAM..."
samtools sort -@ "$CORES" -o "$BAM_FINAL" "$BAM_UNSORTED"
samtools index "$BAM_FINAL"
rm -f "$BAM_UNSORTED"
log "  → $BAM_FINAL (.bai index created)"

# ---------------------------------------------------------------------------
# Step 6: Index .pairs with pairix (required by cooler cload pairix)
# ---------------------------------------------------------------------------
log "Step 6/7: Indexing .pairs with pairix..."
pairix -f "$PAIRS_GZ"
log "  → ${PAIRS_GZ%.gz}.px2"

# ---------------------------------------------------------------------------
# Step 7: Generate .cool and .mcool contact matrices
#
# cooler cload pairix: loads pairs at the specified bin size
# cooler zoomify --balance: generates multi-resolution mcool with ICE
#   normalisation at each zoom level
#
# Key parameter choices:
#   BIN_SIZE=1000  → 1 kb for loop resolution; use 5000 for TAD-only work
#   --balance      → ICE normalisation; always include for downstream TAD calls
#   zoomify creates resolutions: 1kb, 2kb, 5kb, 10kb, 25kb, 50kb, 100kb, ...
#   This library's loader reads any of these via get_cooler(filename, resolution)
# ---------------------------------------------------------------------------
log "Step 7/7: Generating contact matrices..."
cooler cload pairix \
    -p "$CORES" \
    "${CHROMSIZES}:${BIN_SIZE}" \
    "$PAIRS_GZ" \
    "$COOL_OUT"

log "  → $COOL_OUT (single resolution)"

cooler zoomify \
    --balance \
    --nproc "$CORES" \
    -o "$MCOOL_OUT" \
    "$COOL_OUT"

rm -f "$COOL_OUT"  # .mcool supersedes it
log "  → $MCOOL_OUT (multi-resolution, ICE normalised)"

# ---------------------------------------------------------------------------
# QC summary
# ---------------------------------------------------------------------------
log "=== QC Summary ==="
if [[ -f "$DEDUP_STATS" ]]; then
    log "Duplication stats:"
    grep -E '^(total|no_dup|duplicate|cis_1kb)' "$DEDUP_STATS" | \
        awk '{printf "  %-30s %s\n", $1, $2}' || true
fi

PAIRS_COUNT=$(bgzip -d -c "$PAIRS_GZ" 2>/dev/null | grep -v '^#' | wc -l || echo "?")
log "Total no-dup pairs: ${PAIRS_COUNT}"
log ""
log "Thresholds to verify:"
log "  no-dup cis pairs ≥ 1 kb  > 40%  of total no-dup pairs"
log "  total no-dup pairs        > 125M for TAD/loop analysis"
log ""
log "Output .mcool: $MCOOL_OUT"
log "To use with this library:"
log "  from src.loaders import get_cooler"
log "  clr = get_cooler('${MCOOL_OUT}', resolution=5000)"

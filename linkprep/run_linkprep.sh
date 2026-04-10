#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CONFIG_PATH="${1:-$ROOT_DIR/linkprep/config.sh}"
PREPROCESS_SCRIPT="$ROOT_DIR/scripts/preprocess_microc.sh"
QC_SUMMARY="$ROOT_DIR/linkprep/qc_summary.py"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

require_file() {
    local path="$1"
    [[ -f "$path" ]] || {
        echo "ERROR: required file not found: $path" >&2
        exit 1
    }
}

require_tool() {
    local tool="$1"
    command -v "$tool" >/dev/null 2>&1 || {
        echo "ERROR: required tool not found in PATH: $tool" >&2
        exit 1
    }
}

if [[ ! -f "$CONFIG_PATH" ]]; then
    cat >&2 <<EOF
ERROR: config file not found: $CONFIG_PATH

Create it from the example:
  cp "$ROOT_DIR/linkprep/config.example.sh" "$ROOT_DIR/linkprep/config.sh"
EOF
    exit 1
fi

# shellcheck source=/dev/null
source "$CONFIG_PATH"

: "${SAMPLE:?SAMPLE must be set in config}"
: "${R1_FASTQ:?R1_FASTQ must be set in config}"
: "${R2_FASTQ:?R2_FASTQ must be set in config}"
: "${GENOME_FA:?GENOME_FA must be set in config}"
: "${CHROMSIZES:?CHROMSIZES must be set in config}"
: "${OUTPUT_DIR:?OUTPUT_DIR must be set in config}"

require_file "$PREPROCESS_SCRIPT"
require_file "$QC_SUMMARY"
require_file "$R1_FASTQ"
require_file "$R2_FASTQ"
require_file "$GENOME_FA"

if [[ ! -f "$CHROMSIZES" ]]; then
    if [[ -f "${GENOME_FA}.fai" ]]; then
        log "Generating chromosome sizes from ${GENOME_FA}.fai"
        cut -f1,2 "${GENOME_FA}.fai" > "$CHROMSIZES"
    else
        echo "ERROR: missing $CHROMSIZES and ${GENOME_FA}.fai" >&2
        exit 1
    fi
fi

for tool in bwa samtools pairtools pairix cooler bgzip python; do
    require_tool "$tool"
done

mkdir -p "$OUTPUT_DIR"

log "Running LinkPrep preprocessing for sample: $SAMPLE"
log "Output directory: $OUTPUT_DIR"
bash "$PREPROCESS_SCRIPT"

STATS_PATH="$OUTPUT_DIR/${SAMPLE}.dedup.stats"
PAIRS_GZ="$OUTPUT_DIR/${SAMPLE}.mapped.pairs.gz"
MCOOL_PATH="$OUTPUT_DIR/${SAMPLE}.mcool"
QC_MD="$OUTPUT_DIR/${SAMPLE}.qc_summary.md"
QC_JSON="$OUTPUT_DIR/${SAMPLE}.qc_summary.json"

if [[ -f "$STATS_PATH" ]]; then
    log "Writing QC summaries"
    python "$QC_SUMMARY" \
        --stats "$STATS_PATH" \
        --pairs "$PAIRS_GZ" \
        --mcool "$MCOOL_PATH" \
        --sample "$SAMPLE" \
        --markdown-out "$QC_MD" \
        --json-out "$QC_JSON"
else
    log "Skipping QC summary because stats file was not produced"
fi

log "LinkPrep run complete"
log "Primary output: $MCOOL_PATH"
if [[ -f "$QC_MD" ]]; then
    log "QC markdown: $QC_MD"
fi
if [[ -f "$QC_JSON" ]]; then
    log "QC JSON: $QC_JSON"
fi

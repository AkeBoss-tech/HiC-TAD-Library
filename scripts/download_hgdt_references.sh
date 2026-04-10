#!/usr/bin/env bash
# Download hg38 reference files for HG-DT (Phase 2).
# Usage: from repo root,  bash scripts/download_hgdt_references.sh
#
# Environment:
#   HGDT_SKIP_HG38_FA=1   Skip the large ~1 GB hg38.fa.gz download (GENCODE + cCRE still fetch).

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${HGDT_REF_DIR:-$ROOT/data/references}"
mkdir -p "$DEST"

echo "HG-DT references → $DEST"

curl_retry=(curl -fsSL --retry 3 --connect-timeout 30)

# --- GENCODE 49 (annotation GTF) ---
GTF_GZ_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
GTF_GZ="$DEST/gencode.v49.annotation.gtf.gz"
GTF_OUT="$DEST/gencode.v49.annotation.gtf"
if [[ -f "$GTF_OUT" && -s "$GTF_OUT" ]]; then
  echo "✓ GENCODE GTF already present ($(wc -c < "$GTF_OUT") bytes)"
else
  echo "Downloading GENCODE GTF…"
  "${curl_retry[@]}" -o "$GTF_GZ" "$GTF_GZ_URL"
  gunzip -f "$GTF_GZ"
  echo "✓ $GTF_OUT"
fi

# --- SCREEN Registry v4: candidate enhancers (ELS) — single BED, ~94 MB ---
# (Combined GRCh38-cCREs.bed.gz URL from older docs often 404s; ELS is a standard subset.)
ELS_URL="https://downloads.wenglab.org/Registry-V4/GRCh38-cCREs.ELS.bed"
ELS_OUT="$DEST/GRCh38-cCREs.ELS.bed"
if [[ -f "$ELS_OUT" && -s "$ELS_OUT" ]]; then
  echo "✓ SCREEN ELS BED already present ($(wc -c < "$ELS_OUT") bytes)"
else
  echo "Downloading SCREEN GRCh38 ELS cCRE BED (large)…"
  "${curl_retry[@]}" -o "$ELS_OUT" "$ELS_URL"
  echo "✓ $ELS_OUT"
fi

# --- UCSC hg38 primary assembly (very large) ---
HG38_URL="https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
HG38_GZ="$DEST/hg38.fa.gz"
HG38_FA="$DEST/hg38.fa"

if [[ "${HGDT_SKIP_HG38_FA:-}" == "1" ]]; then
  echo "Skipping hg38 FASTA (HGDT_SKIP_HG38_FA=1). Set HGDT_HG38_FASTA manually or re-run without skip."
  exit 0
fi

if [[ -f "$HG38_FA" && -s "$HG38_FA" ]]; then
  echo "✓ hg38.fa already present ($(wc -c < "$HG38_FA") bytes)"
  exit 0
fi

echo "Downloading UCSC hg38.fa.gz (≈1 GB; resumable)…"
# -C - resume partial downloads
curl -fSL --retry 3 --connect-timeout 30 -C - -o "$HG38_GZ" "$HG38_URL" || {
  echo "hg38 download failed or interrupted. Re-run this script to resume."
  exit 1
}
echo "Decompressing hg38.fa.gz…"
gunzip -f "$HG38_GZ"
echo "✓ $HG38_FA"
echo ""
echo "Done. Point HG-DT at these files (defaults match):"
echo "  export HGDT_REF_DIR=$DEST"
echo "  # or: HGDT_HG38_FASTA  HGDT_GENCODE_GTF  HGDT_CCRE_BED"

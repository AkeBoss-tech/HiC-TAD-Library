#!/usr/bin/env bash
#SBATCH --job-name=hichip_download
#SBATCH --partition=main
#SBATCH --account=__LAB_ACCOUNT__
#SBATCH --output=__LOGS_DIR__/01_download_%j.out
#SBATCH --error=__LOGS_DIR__/01_download_%j.err
#SBATCH --time=06:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=__EMAIL__
# =============================================================================
# 01_download_data.slurm
#
# Downloads all Dovetail public HiChIP FASTQ datasets from S3.
# Runs on a single node (bandwidth-limited, no parallelism needed).
#
# Submit:
#   sbatch linkprep/amarel/01_download_data.slurm
# =============================================================================

set -euo pipefail
source "$(dirname "$0")/config.sh"

echo "=== Downloading HiChIP datasets ==="
echo "Destination: ${DATA_DIR}"
date

# Update SLURM header placeholders at runtime (the __VARS__ above are
# for documentation — the actual values come from config.sh sourced above)

mkdir -p "${DATA_DIR}"

# Download each file; skip if already present and complete
download_file() {
    local url="$1"
    local dest="${DATA_DIR}/$(basename "$1")"
    if [[ -f "${dest}" ]]; then
        echo "  SKIP (exists): $(basename "$1")"
        return 0
    fi
    echo "  Downloading: $(basename "$1") …"
    wget -q --show-progress -O "${dest}.tmp" "${url}" && mv "${dest}.tmp" "${dest}"
    echo "  Done: $(basename "$1")  ($(du -sh "$dest" | cut -f1))"
}

for url in "${S3_URLS[@]}"; do
    download_file "${url}"
done

echo
echo "=== Download complete ==="
echo "Files in ${DATA_DIR}:"
ls -lh "${DATA_DIR}"
date

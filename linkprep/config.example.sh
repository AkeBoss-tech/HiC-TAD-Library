#!/usr/bin/env bash
# Example configuration for LinkPrep preprocessing.
# Copy to linkprep/config.sh and edit the values for your dataset.

export SAMPLE="linkprep_sample"

# Paired FASTQ inputs
export R1_FASTQ="data/raw/sample_R1.fastq.gz"
export R2_FASTQ="data/raw/sample_R2.fastq.gz"

# Reference assets
export GENOME_FA="data/raw/mm10.fa"
export CHROMSIZES="data/raw/mm10.chromsizes"

# Output settings
export OUTPUT_DIR="data/processed/linkprep_${SAMPLE}"
export BIN_SIZE="1000"
export CORES="8"
export MIN_MAPQ="40"

# Optional: override TMPDIR if pairtools sort needs a larger scratch disk
# export TMPDIR="/path/to/large/tmp"

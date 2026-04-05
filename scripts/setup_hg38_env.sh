#!/bin/bash
set -e

echo "Setting up HG38 Data Environment..."

# Create target directory
TARGET_DIR="data/reference/hg38"
mkdir -p $TARGET_DIR
cd $TARGET_DIR

echo "Downloading hg38 FASTA..."
rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr*.fa.gz .

echo "Downloading GENCODE v49 GTF..."
wget -qO- https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz > gencode.v49.annotation.gtf.gz || wget https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz

echo "Downloading ENCODE cCREs BED..."
wget -qO- https://screen.wenglab.org/downloads/GRCh38-cCREs.bed.gz > GRCh38-cCREs.bed.gz || wget https://screen.wenglab.org/downloads/GRCh38-cCREs.bed.gz

echo "Done!"

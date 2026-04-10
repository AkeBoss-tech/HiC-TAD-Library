#!/usr/bin/env python3
import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(levelname)s: %(message)s')

def run_cmd(cmd):
    logging.info(f"Running: {cmd}")
    subprocess.run(cmd, shell=True, check=True)

def main():
    target_dir = "data/reference/hg38"
    os.makedirs(target_dir, exist_ok=True)

    # hg38 FASTA
    # using rsync from UCSC
    fasta_cmd = f"rsync -avzP rsync://hgdownload.soe.ucsc.edu/goldenPath/hg38/chromosomes/chr*.fa.gz {target_dir}/"
    run_cmd(fasta_cmd)

    # GENCODE v49 GTF
    gtf_url = "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.annotation.gtf.gz"
    gtf_target = f"{target_dir}/gencode.v49.annotation.gtf.gz"
    gtf_cmd = f"wget -c {gtf_url} -O {gtf_target}"
    run_cmd(gtf_cmd)

    # ENCODE cCREs BED
    bed_url = "https://screen.wenglab.org/downloads/GRCh38-cCREs.bed.gz"
    bed_target = f"{target_dir}/GRCh38-cCREs.bed.gz"
    bed_cmd = f"wget -c {bed_url} -O {bed_target}"
    run_cmd(bed_cmd)

    logging.info("All downloads completed successfully.")

if __name__ == "__main__":
    main()

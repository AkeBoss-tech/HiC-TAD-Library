# HiChIP Pipeline — Rutgers OARC Amarel

End-to-end HiChIP analysis pipeline configured for the Rutgers Amarel cluster
(SLURM scheduler).  Uses the official Dovetail public datasets and scripts from
the [Dovetail-Analysis](https://github.com/dovetail-genomics/Dovetail-Analysis)
submodule at `externals/Dovetail-Analysis/`.

---

## Datasets (Dovetail public, GM12878, hg38)

| Sample | Target | Type | S3 URL |
|--------|--------|------|--------|
| CTCF_2M | CTCF | ~2M reads (test) | `s3://dovetail.pub/HiChIP/fastqs/HiChiP_CTCF_2M_R{1,2}.fastq.gz` |
| CTCF_DS | CTCF | Deep sequencing | `s3://dovetail.pub/HiChIP/fastqs/CTCF-DS_R{1,2}.fastq.gz` |
| H3K27Ac | H3K27Ac | Deep sequencing | `s3://dovetail.pub/HiChIP/fastqs/H3K27Ac_R{1,2}.fastq.gz` |
| H3K4me3 | H3K4me3 | Deep sequencing | `s3://dovetail.pub/HiChIP/fastqs/H3K4me3_R{1,2}.fastq.gz` |

Download with `wget` by replacing `s3://` → `https://s3.amazonaws.com/`
(step 01 does this automatically).

---

## Pipeline overview

```
00_setup_env.sh          ← Run once on login node (conda env + Singularity image)
        │
01_download_data.slurm   ← wget all 8 FASTQ files from S3
        │
02_align.slurm ──────────── BWA-MEM (-5SP -T0), name-sort BAM    [array 0–3]
        │
03_pairtools.slurm ──────── parse → sort → dedup → split          [array 0–3]
        │                   get_qc.py QC summary from stats.txt
        │
   ┌────┴──────────────────┬──────────────────┐
   ▼                       ▼                  ▼
04_contact_matrix.slurm  05_qc.slurm        06_macs2_peaks.slurm
cooler .cool + .mcool    enrichment QC       MACS2 1D peaks
(ICE balanced)           preseq complexity
   │                                          │
   └────────────────────┬─────────────────────┘
                        ▼
               07_fithichip.slurm ──── FitHiChIP loop calling     [array 0–3]
                        │              (Singularity container)
                        ▼
               08_comparative.slurm ── Differential loops (EdgeR)
```

---

## Quick start

### 1. Clone the repo on Amarel

```bash
ssh <netid>@amarel.rutgers.edu
cd /projects/<lab_account>
git clone --recurse-submodules <your_repo_url> hichip_project
cd hichip_project
```

### 2. Edit `config.sh`

```bash
nano linkprep/amarel/config.sh
```

Required changes:
```bash
NETID="your_netid"
LAB_ACCOUNT="your_lab_account"   # check with: sacctmgr show user $USER
```

### 3. One-time setup (login node)

```bash
bash linkprep/amarel/00_setup_env.sh
```

Downloads hg38 (~3 GB), builds BWA index, installs conda env, pulls
FitHiChIP Singularity image, downloads ENCODE peak files.

### 4. Submit the pipeline

```bash
# Full pipeline
bash linkprep/amarel/submit_pipeline.sh

# Dry run first (see what would be submitted)
bash linkprep/amarel/submit_pipeline.sh --dry-run

# Skip download if FASTQs already present
bash linkprep/amarel/submit_pipeline.sh --skip-download
```

### 5. Comparative analysis (after FitHiChIP completes)

```bash
# Edit 08_comparative.slurm: set CONDITION_A and CONDITION_B
sbatch --dependency=afterok:<fithichip_job_id> \
    linkprep/amarel/08_comparative.slurm
```

---

## Amarel tips

### Check your account
```bash
sacctmgr show user $USER withassoc format=User,Account,Partition
```

### Monitor jobs
```bash
squeue -u $USER
sacct -j <jobid> --format=JobID,State,Elapsed,MaxRSS,CPUTime
```

### Storage quotas
```bash
df -h /scratch/$USER     # temporary (not backed up)
df -h /projects/$GROUP   # project (backed up, quota applies)
```

### Interactive session for debugging
```bash
srun --partition=main --cpus-per-task=4 --mem=16G --time=02:00:00 --pty bash
```

---

## QC thresholds (Dovetail)

| Metric | Shallow (20M pairs) | Deep (100–200M pairs) |
|--------|--------------------|-----------------------|
| No-dup read pair rate | > 75% | > 50% |
| No-dup cis ≥ 1 kb | > 20% | > 20% |
| Reads in 1 kb peak window | > 2% | > 2% |
| Merged loops (2.5 kb) | — | > 500 |

---

## Key outputs

| File | Description |
|------|-------------|
| `bam/<sample>.mapped.PT.bam` | Deduplicated, coordinate-sorted BAM |
| `pairs/<sample>.dedup.pairs.gz` | Deduplicated 4DN pairs file |
| `matrices/<sample>.mcool` | Multi-resolution ICE-balanced contact matrix |
| `qc/<sample>.qc_report.txt` | Proximity-ligation QC (get_qc.py) |
| `qc/<sample>.enrichment_qc.txt` | ChIP enrichment statistics |
| `qc/<sample>.enrichment.png` | ChIP enrichment plot |
| `peaks/<sample>.macs2_peaks.bed` | MACS2 1D peak calls |
| `loops/<sample>/*_MergeNearContacts.bed` | FitHiChIP significant loops |
| `loops/differential_*/Loops_EdgeR_Default_SIG.bed` | Differential loops |

---

## References

- Dovetail Analysis docs: https://dovetail-analysis.readthedocs.io/en/latest/hichip/
- FitHiChIP: https://github.com/ay434/FitHiChIP
- pairtools: https://pairtools.readthedocs.io/
- cooler: https://cooler.readthedocs.io/

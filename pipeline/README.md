# Drug Discovery Pipeline

In-silico 4-stage pipeline: **chromatin locus → protein sequence → 3D structure → small-molecule leads**

Built on top of the Unc5b synthesis experiments in `synthesis/`.

## Target

**UNC5B** (Unc-5 Netrin Receptor B) — UniProt Q8K1S3 (mouse)
Death domain, residues 865–943 (~79 aa)
Druggable target in cancer angiogenesis and vascular biology

## Stages

| # | Stage | Tool | API | Est. time |
|---|-------|------|-----|-----------|
| 1 | Locus context + UniProt sequence fetch | UniProt REST | None (public) | ~5 s |
| 2 | ESM-2 protein embeddings | NVIDIA NIM ESM2-30M | `NVIDIA_API_KEY` | ~20 s |
| 3 | ESMFold structure prediction | NVIDIA NIM ESMFold | `NVIDIA_API_KEY` | ~60 s |
| 4 | MolMIM generation + DiffDock docking | NVIDIA NIM | `NVIDIA_API_KEY` | ~10 min |

**Total:** 23 NVIDIA NIM API calls · ~12 min end-to-end

## Quick Start

```bash
# 1. Add your NVIDIA NIM key to .env
echo "NVIDIA_API_KEY=nvapi-..." >> .env

# 2. Run all stages
conda activate hic-analysis
make run-pipeline

# 3. Generate HTML report
make pipeline-report

# 4. Launch Streamlit UI
make pipeline-ui
```

## Individual stages

```bash
make run-pipeline-stage1   # genomics + sequence (no NVIDIA key)
make run-pipeline-stage2   # ESM-2 embeddings
make run-pipeline-stage3   # ESMFold structure
make run-pipeline-stage4   # MolMIM + DiffDock  (~10 min)
```

Force re-run (ignore caches):

```bash
make run-pipeline-force
```

## Outputs

| File | Description |
|------|-------------|
| `data/processed/pipeline_target_sequence.fasta` | UNC5B death domain (~90 aa) |
| `data/processed/pipeline_esm2_embeddings.npy` | ESM-2 embeddings, shape (L, 480) |
| `data/processed/pipeline_structure.pdb` | ESMFold predicted PDB |
| `data/processed/pipeline_molecules.json` | Ranked molecules (rank, SMILES, scores) |
| `media/pipeline_stage1_expression.png` | CTCF signal + gene body locus plot |
| `media/pipeline_stage2_embeddings.png` | Top-20 ESM-2 embedding dimensions |
| `media/pipeline_stage3_structure.png` | 3D Cα backbone trace |
| `media/pipeline_stage4_molecules.png` | DiffDock confidence bar chart (top-10) |
| `pipeline_report.html` | Full HTML report (indigo theme) |

## API Keys Required

- `ALPHA_GENOME_API_KEY` — already in `.env` from synthesis experiments
  (Stage 1 uses the existing WT cache; no new AlphaGenome call needed)
- `NVIDIA_API_KEY` — new key needed for Stages 2–4
  Obtain at: [build.nvidia.com](https://build.nvidia.com)

## Connection to Synthesis Experiments

Stage 1 reads `data/processed/synthesis_unc5b_wt_cache.npz` (created by
`make run-synthesis-exp1`). The CTCF signal context plot shows the same locus
probed in the synthesis experiments, linking chromatin disruption to UNC5B
de-repression and motivating the drug target selection.

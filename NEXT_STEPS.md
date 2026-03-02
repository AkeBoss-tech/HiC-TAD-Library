# Next Steps

## Where We Are

The pipeline has two working layers:

1. **Structural nomination** — Micro-C contact enrichment + TAD membership + boundary exclusion nominates candidate enhancer bins.
2. **CNN sequence scoring** — A dilated CNN (3 × Conv1d with dilation 1→2→4) scores each candidate's 500 bp sequence. Trained on 224 bins using contact enrichment as a proxy for accessibility (loss 0.056→0.004, 30 epochs). Model weights at `data/processed/enhancer_cnn.pt`.

Current figures in `media/` show bars coloured by CNN score (green = high, red = low) with gold ★ stars on called enhancers (score ≥ 0.5).

---

## Step 1 — Retrain CNN on Real ATAC-seq Labels

**Why:** The current model was trained with contact enrichment as a proxy for chromatin accessibility. Real ATAC-seq signal will sharpen specificity considerably — distinguishing actual open chromatin from other contact-rich elements like promoters and CTCF loop anchors.

**What to do:**
- When Jingyun & Ed's new Dovetail Micro-C experiment is sequenced, request that ATAC-seq be run on the same cell type in parallel.
- Convert ATAC-seq BAM → bigWig (e.g. with `bamCoverage --normalizeUsing CPM`).
- In `train_enhancers.py`, replace the `assign_atac_to_candidates()` fallback with a bigWig reader (e.g. `pyBigWig`) that queries signal at each candidate bin's midpoint.
- Retrain with `epochs=50`, `lr=5e-4`; expect loss to converge more slowly on real signal.

**Files to edit:** `train_enhancers.py` (swap label source), `environment.yml` (add `pyBigWig`).

---

## Step 2 — Fix AlphaGenome Mouse ATAC Ontology

**Why:** During training, `CL:0000100` (motor neuron) was rejected for `MUS_MUSCULUS`. AlphaGenome does support mouse ATAC — the ontology term just needs to be correct.

**What to do:**
- Query the AlphaGenome metadata endpoint to list available mouse ATAC tracks:
  ```python
  from alphagenome.models import dna_client, dna_output
  model = dna_client.create(api_key)
  meta = model.get_metadata(organism=dna_model.Organism.MUS_MUSCULUS,
                             output_type=dna_output.OutputType.ATAC)
  print(meta)  # lists available ontology terms
  ```
- Likely candidates: `EFO:0004038` (mouse ESC, already used for deletion scans), `CL:0000540` (neuron), or a UBERON term.
- Once found, update `ATAC_ONTOLOGY` in `train_enhancers.py` and re-run.

---

## Step 3 — Preprocess New Dovetail Micro-C Data

**Why:** Jingyun & Ed's new experiment will provide higher-resolution contact maps to refine candidate calls.

**What to do:**
```bash
# Edit config at top of script first (R1_FASTQ, R2_FASTQ, GENOME_FA, BIN_SIZE)
bash scripts/preprocess_microc.sh
```

**QC thresholds to verify after running:**
- No-dup cis pairs ≥ 1 kb > **40%** of total pairs (ligation quality)
- Total no-dup pairs > **125 M** (depth for TAD/loop calling)

**Then:** update `FILE_PATH` in `visualize_enhancers.py` and `visualize_tads.py` to point to the new `.mcool`.

---

## Step 4 — Expand Candidate Regions

**Why:** The current pipeline covers Sox11 (2 Mb on chr12) and Mir9-2 (1 Mb on chr13). Expanding to the full TAD span or additional loci will give more training data and broader coverage.

**What to do:**
- Add regions to the `regions` dict in `train_enhancers.py` and `visualize_enhancers.py`.
- Jingyun's confirmed full span for the chr13 insulator TAD: `chr13:81,760,002–85,200,000` (3.4 Mb — split into ≤1 Mb windows for AlphaGenome).
- Consider adding a negative-control region (gene desert or heterochromatic locus) to give the CNN clear negative examples.

---

## Step 5 — Synthetic Enhancer Design Loop

**Why:** This is Kelvin's end goal — using the CNN to design novel sequences with high predicted regulatory activity.

**Two approaches from the papers:**

### A. Gradient-based sequence optimization (Cell Systems 2025 approach)
Start from a random or shuffled sequence, compute CNN score, backpropagate into the one-hot input, update sequence toward higher score. Requires making one-hot differentiable (Gumbel-softmax or straight-through estimator).

```python
# Pseudocode
seq = torch.randn(1, 4, 500, requires_grad=True)  # soft one-hot
optimizer = torch.optim.Adam([seq], lr=1e-2)
for step in range(1000):
    score = model(F.softmax(seq, dim=1))
    loss = -score.mean()
    loss.backward()
    optimizer.step()
designed_seq = seq.detach().argmax(dim=1)  # decode to ACGT
```

### B. Motif insertion / combinatorial design
Take a real candidate bin's sequence, mask 8–16 bp windows, and replace with known TF motifs (Sox11, CTCF, AP-1) from JASPAR. Score each variant with the CNN; keep combinations that increase score.

**Files to create:** `design_enhancers.py` (new script), `src/enhancer_design.py` (new module).

---

## Step 6 — Validate Designed Sequences

**Why:** CNN score is a proxy. Experimental validation closes the loop.

**Options (increasing cost):**
1. **AlphaGenome in silico**: run designed sequences through AlphaGenome accessibility predictions to check predicted ATAC signal.
2. **MPRA (Massively Parallel Reporter Assay)**: synthesize a library of designed + control sequences, clone into a reporter vector, transfect into neural cells, read out by sequencing. ~$5–15k, 6–8 week turnaround.
3. **CRISPRi knockdown**: if a designed enhancer overlaps an endogenous candidate, use CRISPRi to repress it and measure Sox11/Mir9-2 expression by RT-qPCR.

---

## Quick Reference — Running the Current Pipeline

```bash
# Nominate + CNN-score candidates (requires mm10.fa + enhancer_cnn.pt)
python visualize_enhancers.py

# Retrain CNN (edit ATAC_ONTOLOGY or label source first)
python train_enhancers.py

# Preprocess new Dovetail data (edit paths at top of script first)
bash scripts/preprocess_microc.sh

# Run all tests
make test-unit
```

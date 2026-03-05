# Next Steps

## Current Pipeline State (as of 2026-03-05)

The full synthetic enhancer pipeline is running end-to-end across three stages:

| Stage | Status | Output |
|-------|--------|--------|
| 1 — Structural nomination (Hi-C) | Complete | 161 Sox11 + 76 Mir9-2 candidate bins |
| 2 — CNN sequence scoring | Complete | All candidates scored; gold ★ bins ≥ 0.5 |
| 3 — Synthetic design | Complete | 20 gradient designs (0.80+), 10 motif-insertion refinements |

**CNN training (2026-03-05):** Real neural ATAC labels from AlphaGenome (forebrain + midbrain + hindbrain, `UBERON:0001890/0001891/0002028`, `MUS_MUSCULUS`). 224 candidates, loss 0.065→0.004, 30 epochs. Model at `data/processed/enhancer_cnn.pt`.

**Design results:**
- Gradient optimisation: best score 0.8155, mean 0.8071, GC 43–47%
- Motif insertion: Sox11 motif accepted in every seed; AP-1, Sox2, NeuroD1 also recurring. Best final score 0.7758 (starting from real top candidate at 0.764). Total improvement per sequence +0.009–+0.050.

---

## Step 1 — AlphaGenome In-Silico Validation of Designs (do now)

**Why:** CNN score alone is a proxy trained on predicted ATAC. Running each designed sequence through AlphaGenome gives an independent second estimate of accessibility — if a designed sequence scores high by CNN *and* by AlphaGenome, confidence in its real-world activity is substantially higher.

**How:** The `validate_designs_alphagenome()` function in `src/enhancer_design.py` is already implemented. It substitutes a designed 500 bp sequence into the real locus context (~1 Mb) and measures delta-ATAC at the substitution site.

**What to write:** `validate_designs.py` (top-level script):
```python
# Pseudocode outline
from src.enhancer_design import validate_designs_alphagenome
from pyfaidx import Fasta

# 1. Load the top-scoring Sox11 locus sequence (2^20 bp centred on chr12:27,120,000)
genome = Fasta("data/raw/mm10.fa")
locus_seq = str(genome["chr12"][26024288:27072864])   # adjust to 2^20 bp

# 2. Load gradient_designs.tsv — take top 10 sequences
designs = pd.read_csv("data/processed/gradient_designs.tsv", sep="\t")["sequence"].head(10).tolist()

# 3. Run AlphaGenome validation (candidate_start = position of top bin within locus_seq)
results = validate_designs_alphagenome(ag_model, locus_seq, locus_interval,
                                        candidate_start=96000,  # ~27,120,000 - 26,024,288
                                        designs=designs,
                                        ontology_terms=["UBERON:0001890", "UBERON:0001891", "UBERON:0002028"])
results.to_csv("data/processed/design_validation.tsv", sep="\t", index=False)
```

**Expected output:** `data/processed/design_validation.tsv` with `delta_atac_mean` and `delta_atac_max` per design. Designs with positive delta are those that increase predicted accessibility above the real sequence — prime candidates for MPRA.

**Files:** create `validate_designs.py`

---

## Step 2 — Expand to More Genomic Regions

**Why:** 224 training examples from 2 loci is a small CNN training set. More regions give:
1. More training data → better generalisation
2. Negative-control regions → teach CNN to distinguish enhancers from other sequence types
3. Broader coverage → more candidate targets for Kelvin's synthetic design goal

**What to add:**

| Region | Coordinates | Purpose |
|--------|-------------|---------|
| Chr13 full TAD span | `chr13:81,760,002-85,200,000` (3.4 Mb) | Expand Mir9-2 coverage (split into 3 × 1 Mb AlphaGenome windows) |
| Negative control (gene desert) | `chr6:90,000,000-91,000,000` (or similar CTCF-only B-compartment) | Give CNN clear negative examples |
| Promoter-control | Pick 3 known active promoters in Sox11 TAD | Compare CNN sensitivity to promoters vs. distal enhancers |

**Files to edit:** `train_enhancers.py` (add to `regions` dict), `visualize_enhancers.py` (same), `design_enhancers.py` (same).

---

## Step 3 — Preprocess New Dovetail Micro-C Data

**Why:** Jingyun & Ed's new Dovetail LinkPrep experiment will provide higher-resolution Micro-C from the exact cell type of interest. The current pipeline uses a downloaded public dataset (mouse brain Micro-C). New data will:
- Confirm which TAD boundaries are real in the actual experimental cell type
- Provide sharper contact maps at 1 kb resolution (vs current 5 kb)
- Enable direct comparison: WT vs. any deletion experiments

**What to do:**
```bash
# 1. Edit config at top of script (R1_FASTQ, R2_FASTQ, GENOME_FA, BIN_SIZE)
vim scripts/preprocess_microc.sh

# 2. Run full pipeline
bash scripts/preprocess_microc.sh

# 3. QC thresholds to verify
# cis pairs ≥ 1 kb > 40% of total pairs
# total no-dup pairs > 125 M
```

**Then:** update `FILE_PATH` in all visualization and training scripts to point to new `.mcool`.

---

## Step 4 — Retrain CNN on Real ATAC-seq Labels

**Why:** The current CNN uses AlphaGenome-*predicted* ATAC as labels. Real ATAC-seq from the experimental cell type will:
- Match the exact chromatin state being studied
- Provide higher-resolution peak calls (ATAC-seq peaks are ~200 bp; our bins are 500 bp — near-perfect match)
- Remove any AlphaGenome prediction artifacts

**What to do:**
1. Run ATAC-seq in parallel with the new Dovetail Micro-C experiment (same cell type, same passage)
2. Process ATAC-seq: trim (Trimmomatic) → align (bowtie2) → remove duplicates (picard) → call peaks (MACS2) → generate bigWig (bamCoverage --normalizeUsing CPM)
3. In `train_enhancers.py`, replace `get_atac_for_region()` with a bigWig reader:
```python
import pyBigWig
def get_atac_from_bigwig(bw_path: str, chrom: str, start: int, end: int) -> np.ndarray:
    bw = pyBigWig.open(bw_path)
    signal = bw.stats(chrom, start, end, type="mean", nBins=(end-start)//128)
    return np.array([v or 0.0 for v in signal], dtype=np.float32)
```
4. Retrain with `epochs=50`, `lr=5e-4`, more stringent peak labels (0/1 from peak calls rather than continuous signal)

**Files to edit:** `train_enhancers.py`, `environment.yml` (add `pyBigWig`)

---

## Step 5 — Build the MPRA Library

**Why:** This is the wet-lab validation step — the only way to know if designed sequences actually drive transcription in neural cells.

**Library composition (~200 sequences, 200 bp inserts for Twist synthesis):**

| Category | Count | Source |
|----------|-------|--------|
| Top gradient designs (Sox11 region) | 20 | `gradient_designs.tsv` — score > 0.80 |
| Top gradient designs (Mir9-2 region) | 20 | `gradient_designs.tsv` |
| Motif-insertion refinements | 10 | `motif_designs.tsv` — best 10 seeds |
| Real top candidates | 20 | Highest CNN-score real bins from both regions |
| Real bottom candidates | 20 | Lowest-score bins (negative controls within TADs) |
| AlphaGenome-validated designs | 10 | Top designs by delta-ATAC from `design_validation.tsv` |
| Positive controls | 5 | Known active enhancers (e.g. Sox11 promoter region, published VISTA elements) |
| Negative controls (scrambled) | 20 | Dinucleotide-shuffled versions of top designs (preserves composition but destroys motifs) |
| Negative controls (gene desert) | 10 | Random sequences from B-compartment gene desert |
| Length-matched random | 15 | Randomly sampled mm10 sequences matched for GC content |

**Primer design:** Each insert flanked by universal amplification primers + barcode (random 20-mer). Read out barcode frequency by sequencing after transfection.

**Cell type recommendation:** Primary mouse cortical neurons or NSC (neural stem cells) — matches Sox11/Mir9-2 expression context. Alternatively, Neuro2a (N2a) mouse neuroblastoma cell line for easier transfection.

**Timeline estimate:** Library design (~1 week) → Twist synthesis (~3 weeks) → cloning into reporter vector (~1 week) → transfection + sequencing (~2 weeks) = ~7 weeks total.

---

## Step 6 — Improve the CNN Architecture

**Why:** The current CNN (3 × Conv1d, dilation 1/2/4, 32 filters) is deliberately shallow. Once more training data is available (Step 2 + Step 4), a more expressive model will capture richer regulatory grammar.

**Options (roughly increasing complexity):**

### A. Deeper dilated CNN (next step)
Add more dilation layers and filter widths to cover a wider effective receptive field:
```python
# Add dilation=8 and dilation=16 layers; increase filters from 32 to 64
# Effective receptive field: 1 + 2*(kernel-1)*(1+2+4+8+16) = ~330 bp for kernel=11
```

### B. Transformer attention over motif embeddings
Replace GlobalAvgPool with a multi-head attention layer. This allows the model to capture long-range motif cooperativity (e.g. Sox11 at position 100 cooperating with AP-1 at position 400).

### C. Pre-trained genomic language model (Nucleotide Transformer / DNABERT2)
Fine-tune a pre-trained genomic transformer on our ATAC labels. These models have seen billions of base pairs and have learned general regulatory grammar. Fine-tuning on 224 (or more) examples is feasible with LoRA.

**Practical recommendation:** Start with (A) — doubling dilation layers is a one-line change and unlikely to overfit even on 224 examples. Move to (C) once real ATAC labels (Step 4) provide more signal.

---

## Step 7 — Iterative Design Loop

**Why:** The Cell Systems 2025 paper uses an iterative loop: design → predict → wet-lab validate → retrain on validated sequences → redesign. Each cycle improves the model's understanding of what makes a real enhancer.

**Loop structure:**
```
Train CNN on available labels (Steps 2, 4)
    ↓
Design synthetic sequences (Step 3 / design_enhancers.py)
    ↓
In-silico filter with AlphaGenome (Step 1 of this doc)
    ↓
Wet-lab MPRA validation (Step 5)
    ↓
Add validated sequences (positive + negative) to training set
    ↓
Retrain CNN on augmented dataset → GOTO top
```

After 2–3 cycles, the CNN will have seen both computationally designed and experimentally validated sequences, substantially improving its specificity for real enhancer activity vs. other contact-rich elements.

---

## Quick Reference — Running the Current Pipeline

```bash
# Full pipeline (requires mm10.fa, AlphaGenome API key, trained model)
python train_enhancers.py          # retrain CNN with real neural ATAC labels
python visualize_enhancers.py      # produce CNN-scored candidate figures
python design_enhancers.py         # gradient + motif-insertion design
python validate_designs.py         # AlphaGenome in-silico validation (to write)

# Preprocessing (when new Dovetail data arrives)
bash scripts/preprocess_microc.sh  # FASTQ → .mcool

# Tests
make test-unit                     # fast unit tests (no data needed)
make test                          # all tests including integration

# Key output files
data/processed/enhancer_cnn.pt          # trained model weights
data/processed/gradient_designs.tsv    # 20 gradient-optimised sequences + scores
data/processed/motif_designs.tsv       # 10 motif-insertion sequences + insertion history
```

---

## AlphaGenome Mouse Track Reference

Confirmed available tracks for `MUS_MUSCULUS` (queried 2026-03-05):

| Track | Ontology | Type | Notes |
|-------|----------|------|-------|
| Forebrain ATAC | `UBERON:0001890` | tissue | **Primary label source** |
| Midbrain ATAC | `UBERON:0001891` | tissue | Averaged with forebrain |
| Hindbrain ATAC | `UBERON:0002028` | tissue | Averaged with forebrain |
| Neural tube ATAC | `UBERON:0001049` | tissue | Alternative — more embryonic |
| Mouse ESC | `EFO:0004038` | cell line | Used in deletion scans; not neural |

The three neural tissue tracks (`UBERON:0001890/0001891/0002028`) are averaged to create a single pan-neural ATAC signal. This is the current CNN training label source.

---

## Key Biological Interpretations

**Sox11 region (chr12:26–28 Mb):**
- 161 candidate bins across 4+ TADs
- Top candidate: chr12:27,120,000–27,125,000 (CNN score 0.764, GC 42%)
- High-scoring bins cluster in the left TAD (~26–26.75 Mb) — the distal regulatory landscape of Sox11
- Sox11 motifs dominate gradient designs and motif insertions — consistent with Sox11 being its own regulator (autoregulation)

**Mir9-2 region (chr13:83.5–84.5 Mb):**
- 76 candidate bins across 2 TADs separated by Jingyun's confirmed insulator
- Top candidate: chr13:84,025,000–84,030,000 (CNN score 0.677, GC 37%)
- The right-side cluster (~84.25–84.4 Mb) shows convergent structural + sequence signal — highest priority for wet-lab testing
- Sox11 + Sox2 motifs dominate insertions — consistent with Mir9-2 being regulated by neural SOX factors

**Gradient design interpretation:**
- All 20 designs score ≥ 0.80 vs. best real candidate at 0.764 — the model has found a region of sequence space predicted to be more accessible than endogenous regulatory elements
- This is expected: gradient optimisation maximises the training objective, so designs should score near the model's saturation point
- The real test is whether AlphaGenome (independent model) and MPRA (wet lab) agree

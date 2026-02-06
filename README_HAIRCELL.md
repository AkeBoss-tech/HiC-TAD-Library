# Hair Cell Differentiation: Hi-C + AlphaGenome Integration

## Biological Question

**How do inner hair cells (IHC) and outer hair cells (OHC) in the mammalian cochlea acquire their distinct cellular identities from the same progenitor cells?**

Both cell types originate from common sensory progenitors that express *ATOH1* (a master regulator of hair cell fate). Yet they develop radically different molecular profiles and functions:

- **Inner Hair Cells (IHC)**: ~3,500 cells in the human cochlea. Primary auditory sensory receptors. Express *SLC17A8* (VGLUT3) for glutamate release to spiral ganglion neurons. Encode frequency and intensity information for sound perception.

- **Outer Hair Cells (OHC)**: ~12,000 cells in three rows. Cochlear amplifiers. Express *SLC26A5* (Prestin), a motor protein that changes cell length in response to voltage. Provide 40-60 dB of mechanical amplification to sharpen frequency tuning.

The central mystery: **what are the "instructions" in the DNA that cause one cell to activate *Prestin* and become an OHC, while its neighbor does not?**

---

## Hypothesis: 3D Genome Architecture Controls Lineage Specification

Gene regulation is not just about transcription factors binding to DNA — it requires that **enhancers physically contact their target gene promoters** through 3D chromatin looping. These loops are constrained by **Topologically Associating Domains (TADs)**, which are structural "neighborhoods" bounded by CTCF/cohesin proteins.

We hypothesize that:

1. **Some TAD boundaries are "hard-coded" in the DNA sequence** (via CTCF motifs) and are conserved between species → these define the *hardware* of genome organization.

2. **Other boundaries are cell-type-specific** and depend on epigenetic modifications → these define the *software* that rewires connectivity during differentiation.

3. **Mutations that disrupt sequence-encoded boundaries** can cause genes to "escape" into neighboring regulatory domains, leading to ectopic expression and hearing loss.

By comparing **observed 3D structure** (from mouse Micro-C Hi-C data) with **predicted structure** (from AlphaGenome's DNA sequence-based model), we can identify which boundaries are intrinsic to the sequence vs. which require cellular context — and then test which sequence changes would abolish hair cell identity.

---

## Overview of the Analysis Pipeline

This workflow integrates three complementary data sources:

| Source | Type | Resolution | Provides |
|--------|------|------------|----------|
| **Mouse Micro-C** (observed) | Experimental chromatin contacts | 5 kb bins | Ground-truth TAD boundaries in *in vivo* tissue |
| **AlphaGenome** (predicted) | AI model trained on 5,930 human epigenomic tracks | 128 bp resolution | Cell-type-specific predictions of contacts, TF binding, accessibility, expression |
| **LiftOver** (mm10→hg38) | Sequence homology mapping | Single-nucleotide | Translation layer to compare orthologous loci |

The pipeline consists of **five steps** that progressively narrow from genome-scale structure to single-nucleotide functional variants:

```
Step 1: Boundary Detection (Mouse Hi-C)
   ↓
Step 2: Cell-Type Predictions (AlphaGenome: IHC vs OHC)
   ↓
Step 3: Boundary Validation (Sequence-Encoded vs Epigenetic)
   ↓
Step 4: Virtual 4C Loop Analysis (Enhancer-Promoter Contacts)
   ↓
Step 5: In Silico Mutagenesis (Variant Prioritization)
```

---

## Step-by-Step Biological Logic

### **Step 1: Detect TAD Boundaries in Mouse Micro-C and Lift to Human**

**Script**: [`haircell_workflow.py::step1_detect_and_liftover()`](haircell_workflow.py)
**Calls**: [`src/tad_boundaries.py`](src/tad_boundaries.py), [`src/liftover_utils.py`](src/liftover_utils.py)

#### What it does:
1. Load mouse Micro-C contact matrices at 5kb resolution around hair cell gene loci (*Atoh1*, *Slc26a5*, *Pou4f3*, *Gfi1*)
2. Run **multi-scale insulation scoring** (diamond-window approach) at 8 logarithmically-spaced window sizes (25kb–500kb)
3. Detect **insulation minima** using topographic prominence (scipy peak-finding)
4. Classify boundaries as:
   - **Strong** (prominence ≥ 0.5) → likely CTCF-anchored, stable across cell types
   - **Weak** (0.2–0.5) → may be transient or sub-TAD boundaries
   - **Sub-threshold** (< 0.2) → noise
5. Convert mouse mm10 coordinates to human hg38 using `pyliftover`

#### Biological interpretation:
- **Strong boundaries** that survive liftover are candidates for sequence-encoded structural anchors — these should be conserved between mouse and human, and AlphaGenome should predict them from sequence alone.
- **Weak boundaries** may be tissue-specific or developmental-stage-specific.
- Boundaries that fail liftover may be species-specific regulatory innovations.

#### Output:
- `media/haircell/{Gene}_lifted_boundaries.csv` — table with source/target coordinates and liftover success flags

---

### **Step 2: Predict IHC and OHC Regulatory Profiles with AlphaGenome**

**Script**: [`haircell_workflow.py::step2_predict_haircell_profiles()`](haircell_workflow.py)
**Calls**: AlphaGenome API (`dna_client.predict_interval`)

#### What it does:
1. For each human gene locus (1 Mb window), request AlphaGenome predictions using cell-type ontology terms:
   - **IHC**: `CL:0000202` (inner hair cell)
   - **OHC**: `CL:0000601` (outer hair cell)
   - **Fallback**: `UBERON:0001846` (internal ear, if specific cell types unavailable)

2. Request four modalities:
   - `CONTACT_MAPS` — predicted TAD structure (128 bp × 128 bp matrix)
   - `RNA_SEQ` — predicted gene expression levels (strand-specific)
   - `ATAC` — predicted chromatin accessibility
   - `CHIP_TF` — predicted transcription factor binding (focus on **CTCF** for boundaries)

3. Compute **differential contact map**: IHC − OHC
   - Positive values (blue) = regions with stronger contacts in IHC
   - Negative values (red) = regions with stronger contacts in OHC

4. Extract and compare **mean CTCF signal** between IHC and OHC

#### Biological interpretation:
- If a TAD boundary **exists in both IHC and OHC** (same CTCF peaks, same insulation dip), it is likely a constitutive structural boundary.
- If a boundary **appears only in OHC** or has differential CTCF strength, it may regulate cell-type-specific enhancer-promoter wiring.
- Differential contacts reveal **domains that physically reorganize** during differentiation — these may contain lineage-specific enhancers that gain access to *Prestin* only in OHCs.

#### Output:
- `media/haircell/{Gene}_diff_contacts.png` — heatmap of IHC − OHC contact differences
- `media/haircell/{Gene}_ctcf_comparison.png` — IHC vs OHC CTCF signal tracks

---

### **Step 3: Validate Which Boundaries Are Sequence-Encoded**

**Script**: [`haircell_workflow.py::step3_validate_boundaries()`](haircell_workflow.py)
**Calls**: [`src/haircell_integration.py`](src/haircell_integration.py)

#### What it does:
1. Compute **insulation score** from AlphaGenome's predicted contact map using the same diamond-window algorithm as Step 1
2. Detect **predicted boundaries** (insulation minima with prominence scoring)
3. For each lifted mouse boundary, find the nearest predicted boundary within ±10 kb
4. Classify concordance:
   - **Concordant** (distance ≤ tolerance) → **"sequence-encoded"** — the boundary exists because of DNA sequence features (CTCF motifs) that AlphaGenome recognized
   - **Discordant** (no match) → **"epigenetically regulated"** — the boundary exists in real cells but AlphaGenome missed it, implying it requires cell-state-specific factors not encoded in the DNA sequence

#### Biological interpretation:
This is the **key validation step** that distinguishes:

| Boundary Type | In Mouse Hi-C? | In AlphaGenome? | Mechanism |
|---------------|----------------|-----------------|-----------|
| **Sequence-encoded** | ✅ | ✅ | CTCF motifs create intrinsic insulators; conserved across species and cell types |
| **Epigenetically regulated** | ✅ | ❌ | Requires pioneer factors, DNA methylation, or histone marks to establish; may be hair-cell-specific |
| **False positive** | ❌ | ✅ | Predicted from sequence but not biologically active (model artifact) |
| **Species-specific** | ✅ | N/A (liftover failed) | Regulatory innovation unique to mouse or human |

**Clinical relevance**: Variants in sequence-encoded boundaries are more likely to cause **dominant genetic hearing loss** (e.g., DFNA mutations) because disrupting a CTCF motif will abolish the boundary in all cell types. Variants in epigenetically regulated boundaries may cause subtler, incomplete, or recessive phenotypes.

#### Output:
- `media/haircell/{Gene}_concordance.csv` — per-boundary match table
- `media/haircell/{Gene}_boundary_validation.png` — insulation track with color-coded boundaries (green = conserved, red = epigenetic)

---

### **Step 4: Virtual 4C — Map Enhancer-Promoter Loops**

**Script**: [`haircell_workflow.py::step4_virtual_4c()`](haircell_workflow.py)
**Calls**: [`src/haircell_integration.py::extract_virtual_4c()`](src/haircell_integration.py), AlphaGenome `to_virtual_4c()`

#### What it does:
1. Set the **viewpoint** at the gene promoter (e.g., *SLC26A5* TSS at chr7:107,301,550)
2. Extract one row of the contact matrix = **virtual 4C profile** (all interactions with the promoter)
3. Compute separately for IHC and OHC predicted contacts
4. Identify **cell-type-specific loops**:
   - Peaks in `IHC − OHC` profile (>95th percentile) = **OHC-specific enhancers** that contact the promoter only in outer hair cells
   - Peaks in `OHC − IHC` profile = **IHC-specific enhancers**

5. Generate AlphaGenome's native virtual 4C visualization with gene transcript annotations

#### Biological interpretation:
This step answers: **"Which regulatory elements physically contact this gene in each cell type?"**

Example for *SLC26A5* (Prestin):
- In **OHC**, we expect to see a strong loop connecting a distal enhancer (e.g., +50 kb downstream) to the *Prestin* promoter.
- In **IHC**, that loop should be absent or weaker, explaining why *Prestin* is not expressed.
- The enhancer region should have:
  - High ATAC signal (open chromatin) in OHC only
  - OHC-specific TF binding (e.g., GATA3, which is enriched in OHCs)
  - Fall within the same TAD as *Prestin* (no boundary blocking it)

**Clinical application**: Variants in OHC-specific enhancers can cause **non-syndromic hearing loss** without affecting IHC function — patients have impaired amplification (threshold elevation) but intact primary transduction.

#### Output:
- `media/haircell/{Gene}_v4c_comparison.png` — overlay of IHC and OHC virtual 4C with difference plot
- `media/haircell/{Gene}_v4c_alphagenome.png` — AlphaGenome native plot with transcript annotations
- Console output listing genomic coordinates of top cell-type-enriched interaction bins

---

### **Step 5: In Silico Mutagenesis — Prioritize Pathogenic Variants**

**Script**: [`haircell_workflow.py::step5_ism_boundaries()`](haircell_workflow.py)
**Calls**: AlphaGenome `predict_variant()`, [`src/haircell_integration.py`](src/haircell_integration.py)

#### What it does:
1. Select **sequence-encoded boundaries** from Step 3 (concordant matches)
2. For each boundary, locate the **CTCF binding peak** within ±5 kb
3. At the CTCF peak position, generate all 12 possible single-nucleotide variants:
   ```
   A→C, A→G, A→T
   C→A, C→G, C→T
   G→A, G→C, G→T
   T→A, T→C, T→G
   ```
4. For each variant, call `predict_variant()` to get:
   - **Reference contact map** (wild-type sequence)
   - **Alternate contact map** (with the mutation)
5. Compute insulation score for both maps
6. Calculate **TAD Dissolution Score**:
   ```
   score = mean(|insulation_alt − insulation_ref|)
   ```
   - High score = variant strongly disrupts boundary insulation
   - Low score = variant has minimal structural effect

7. Rank variants by score and identify the top 5 most disruptive

#### Biological interpretation:
This is a **saturating mutagenesis** of the CTCF motif to identify which bases are most critical for boundary function.

**Example output** (hypothetical):
```
Top variants by TAD dissolution:
  G>A at 107,850,234  score = 0.428   ← destroys core CTCF motif
  C>T at 107,850,237  score = 0.391
  G>C at 107,850,235  score = 0.312
  T>G at 107,850,240  score = 0.089   ← flanking base, minor effect
  A>C at 107,850,245  score = 0.012   ← no effect
```

These predictions can be cross-referenced with:
- **ClinVar / gnomAD**: Are these positions already linked to hearing loss?
- **CTCF PWM (position weight matrix)**: Do high-scoring variants match critical positions in the CTCF binding consensus?
- **Conservation (PhyloP/GERP)**: Are these bases constrained across vertebrates?

**Mechanistic hypothesis**:
If a variant disrupts a boundary upstream of *SLC26A5*, the gene may:
1. Lose access to its OHC-specific enhancer (if the enhancer is now in a different TAD) → **loss of Prestin expression** → OHC dysfunction
2. Gain access to a silencer or inappropriate enhancer from an adjacent TAD → **ectopic expression** or **dysregulation**

#### Output:
- `media/haircell/{Gene}_ism_scores.csv` — ranked table of all tested variants
- `media/haircell/{Gene}_ism_scores.png` — bar chart with top 5 variants annotated
- Console output with ref>alt notation and scores

---

## How the Scripts Work Together

### File Organization

```
HiC-TAD-Library/
├── haircell_workflow.py              ← MAIN SCRIPT (orchestrates Steps 1–5)
├── src/
│   ├── tad_boundaries.py             ← Hi-C analysis (insulation, prominence, multi-scale)
│   ├── liftover_utils.py             ← Coordinate conversion mm10→hg38
│   ├── haircell_integration.py       ← Integration functions + plotting
│   └── loaders.py                    ← Cooler file loading
├── visualize_boundaries.py           ← Standalone boundary detection visualizations
├── visualize_alphagenome.py          ← Standalone AlphaGenome prediction plots
├── tal1_example_workflow.py          ← Template for variant scoring (TAL1 oncogene)
├── data/raw/mouse_microc.mcool       ← 4.8 GB mouse Micro-C contact matrices
├── media/haircell/                   ← All outputs go here
└── .env                              ← ALPHA_GENOME_API_KEY
```

### Data Flow

```
Mouse Micro-C (mm10)
   └─> [tad_boundaries.py] → Insulation scores, boundary positions
         └─> [liftover_utils.py] → Human hg38 coordinates
               │
Human Gene Loci (hg38)                                AlphaGenome API
   └─> [haircell_workflow.py]                             │
         ├─ Step 2 → predict_interval(IHC/OHC) ←─────────┘
         │            ├─ Contact maps (predicted TADs)
         │            ├─ CTCF binding
         │            ├─ ATAC accessibility
         │            └─ RNA-seq
         │
         ├─ Step 3 → [haircell_integration.py]
         │            ├─ Insulation from predicted contacts
         │            └─ Concordance analysis
         │
         ├─ Step 4 → extract_virtual_4c()
         │            └─ Loop identification
         │
         └─ Step 5 → predict_variant() ←──────────────────┘
                      └─ Variant effect scores
                           └─ [haircell_integration.py]
                                └─ Plot ISM rankings
```

### Module Responsibilities

| Module | Role |
|--------|------|
| **`haircell_workflow.py`** | High-level orchestration; loads data, calls all steps, handles I/O |
| **`tad_boundaries.py`** | Core Hi-C algorithms (directionality index, insulation, prominence scoring, TAD calling) |
| **`liftover_utils.py`** | Genomic coordinate translation using `pyliftover` |
| **`haircell_integration.py`** | Unified analysis layer (works on both observed and predicted data); boundary concordance; virtual 4C extraction; all visualization functions |
| **`loaders.py`** | Wrapper for cooler file access |

---

## Running the Analysis

### Prerequisites

```bash
# Install dependencies
make install

# Required files:
# - data/raw/mouse_microc.mcool (4.8 GB, run notebooks/00_fetch_data.py)
# - .env with ALPHA_GENOME_API_KEY=your_key
```

### Basic Usage

```bash
# Run for one gene
python haircell_workflow.py ATOH1

# Or via make
make run-haircell GENE=SLC26A5

# Run all four genes
python haircell_workflow.py all
```

### Expected Runtime

| Step | Time (per gene) | Notes |
|------|----------------|-------|
| Step 1 | ~30 sec | Local Hi-C analysis |
| Step 2 | ~2–5 min | Two AlphaGenome API calls (1 MB intervals, 4 modalities each) |
| Step 3 | ~10 sec | Local analysis of predicted contacts |
| Step 4 | ~1 min | Virtual 4C + AlphaGenome native plot |
| Step 5 | ~10–30 min | 12 SNVs × N boundaries (API rate-limited) |

**Total**: ~15–40 minutes per gene (dominated by ISM variant scoring)

If Step 5 fails (CONTACT_MAPS not available for `predict_variant`), the script automatically falls back to RNA_SEQ-based scoring.

---

## Interpreting the Results

### Key Output Files

#### 1. Boundary Liftover (`{Gene}_lifted_boundaries.csv`)
Columns:
- `source_chrom`, `source_start`, `source_end` — mouse mm10 coordinates
- `target_chrom`, `target_pos` — human hg38 coordinate
- `boundary_class` — strong/weak
- `prominence` — insulation dip depth
- `liftover_success` — TRUE if mapped

**What to look for**:
- Strong boundaries with high liftover success rate → conserved structural anchors
- Failed liftovers → species-specific regulation

#### 2. Concordance Table (`{Gene}_concordance.csv`)
Columns:
- `observed_bin` — lifted mouse boundary position
- `nearest_predicted_bin` — AlphaGenome predicted boundary position
- `distance_bins` — separation
- `is_concordant` — match within tolerance?
- `interpretation` — "sequence-encoded" or "epigenetically regulated"

**What to look for**:
- High concordance rate (>70%) → locus has strong sequence-encoded structure
- Low concordance → cell-type-specific or developmental regulation dominates

#### 3. ISM Scores (`{Gene}_ism_scores.csv`)
Columns:
- `boundary_pos` — genomic location of boundary
- `variant_pos` — CTCF peak position
- `ref`, `alt` — nucleotide change
- `tad_dissolution_score` — mean insulation change

**What to look for**:
- Variants with score >0.3 → strong boundary disruption
- Check if high-scoring variants fall in conserved CTCF motif positions (JASPAR database)
- Cross-reference with hearing loss mutation databases

### Visual Outputs

#### Differential Contact Map
**Blue regions**: Stronger contacts in IHC
**Red regions**: Stronger contacts in OHC
**Interpretation**: Domains that physically reorganize during differentiation

#### Boundary Validation Plot
**Green dashed lines**: Sequence-encoded boundaries (conserved)
**Red dashed lines**: Epigenetically regulated boundaries
**Interpretation**: Green = targetable by sequence variants, Red = requires epigenetic perturbation

#### Virtual 4C Comparison
**Top panel**: Overlay of IHC (blue) and OHC (orange) contact profiles
**Bottom panel**: Difference track (IHC − OHC)
**Interpretation**: Peaks in difference track = cell-type-specific enhancers

#### ISM Bar Chart
**Red bars**: Top 5 most disruptive variants
**Blue bars**: Minor-effect variants
**Interpretation**: Red variants are candidates for functional validation (CRISPR)

---

## Biological Predictions & Testable Hypotheses

### Prediction 1: Prestin Locus Has an OHC-Specific Enhancer Loop

**Test**: Compare *SLC26A5* virtual 4C output between IHC and OHC. We predict:
- A strong peak ~50–150 kb from the TSS in OHC only
- That region should have high ATAC signal in OHC (from Step 2 output)
- A sequence-encoded boundary should insulate this enhancer from neighboring genes

**Validation experiment**:
- CRISPR-delete the predicted enhancer in mouse cochlear organoids
- Measure *Slc26a5* expression by RNA-FISH or qPCR
- Expected result: Loss of Prestin in OHCs, no effect on IHCs

### Prediction 2: ATOH1 Boundary Mutations Cause Ectopic Gene Activation

**Test**: Look at `ATOH1_ism_scores.csv` for high-scoring variants in boundaries flanking *ATOH1*.

**Validation experiment**:
- Introduce top-scoring variant into human iPSC-derived otic progenitors via base editing
- Differentiate to hair cells
- Perform RNA-seq to detect genes that become mis-regulated
- Expected result: Genes from adjacent TADs gain inappropriate contact with *ATOH1* enhancers

### Prediction 3: Some Boundaries Are Epigenetically Regulated by Pioneer Factors

**Test**: Boundaries marked "epigenetically regulated" in `{Gene}_concordance.csv` should have:
- Low CTCF signal in AlphaGenome predictions
- Enrichment for pioneer TF binding sites (SOX2, PAX2, GFI1)

**Validation experiment**:
- ChIP-seq for candidate pioneer factors in sorted hair cells
- Compare binding at epigenetic vs sequence-encoded boundaries
- Expected result: Epigenetic boundaries recruit non-CTCF architectural proteins

---

## Limitations & Caveats

### 1. AlphaGenome Training Data Bias
- Model is trained on bulk human tissues, not single-cell hair cell data
- IHC/OHC ontology terms may map to broader "inner ear" tracks if specific cell types are missing
- Prediction accuracy for rare cell types is untested

### 2. Cross-Species Comparison
- Mouse-to-human liftover assumes regulatory conservation (not always true)
- Species-specific TAD structure has been documented (e.g., *EPHA4* locus)
- Mouse cochlea has 3 rows of OHCs, humans have variable numbers

### 3. Insulation Score Sensitivity
- Diamond window size (50 kb default) affects boundary detection
- Boundaries at multiple scales may represent nested sub-TADs vs true insulators
- Prominence thresholds (0.2/0.5) are empirical and dataset-dependent

### 4. ISM Incompleteness
- Only tests SNVs at CTCF peaks, not insertions/deletions or structural variants
- Does not account for compound heterozygous effects
- `predict_variant` with CONTACT_MAPS may not be supported; falls back to RNA_SEQ scoring

### 5. No Ground-Truth Hair Cell Hi-C
- Mouse Micro-C is from whole tissue (mixed cell types)
- Ideally, we would have sorted IHC and OHC Hi-C for direct comparison
- Predictions remain to be validated experimentally

---

## Future Extensions

### 1. Integrate Developmental Time Series
- Run workflow on *Atoh1*-null vs wild-type Hi-C to see how master regulators control boundaries
- Predict when during development (E13.5, E16.5, P0, P7) boundaries form

### 2. Allele-Specific Analysis
- For hearing loss patients with *SLC26A5* variants, run ISM on their exact mutation
- Compare predicted TAD dissolution score to audiogram phenotype

### 3. Comparative Genomics
- Liftover to zebrafish, chicken, bat (echolocating species with extreme OHC specialization)
- Identify lineage-specific boundary innovations

### 4. CRISPR Validation
- Use boundary predictions to design CTCF motif disruptions
- Insert reporters (GFP) downstream of predicted cell-type-specific enhancers
- Perform live imaging in cochlear organoids

### 5. Drug Target Identification
- If an epigenetically regulated boundary requires BRD4 (bromodomain protein), test BET inhibitors
- Predict which small molecules could "flip" a boundary to reactivate lost *Prestin* expression

---

## References & Acknowledgments

### Key Methods Papers

**TAD Boundary Detection**:
- Dixon et al. (2012) *Nature* — Original directionality index method
- Crane et al. (2015) *Nature* — Insulation score diamond-window approach
- Rao et al. (2014) *Cell* — Hi-C at kb resolution, loop detection

**AlphaGenome**:
- Kelley et al. (2018) *Genome Research* — Enformer (precursor model)
- Avsec et al. (2021) *Nature Methods* — Basenji architecture
- *AlphaGenome Documentation* — https://www.alphagenomedocs.com/

**Hair Cell Biology**:
- Bermingham et al. (1999) *Science* — *ATOH1* is required for hair cell fate
- Zheng et al. (2000) *Nature* — *Prestin* (SLC26A5) is the OHC motor protein
- Kolla et al. (2020) *eLife* — GFI1 controls IHC vs OHC differentiation

**3D Genome & Development**:
- Bonev et al. (2017) *Cell* — TAD formation during mouse embryogenesis
- Kragesteen et al. (2018) *Nature* — Limb enhancer-promoter specificity via TAD boundaries

### Data Sources

- **Mouse Micro-C**: 4DNucleome Consortium (4DNES14CNC1I)
- **GENCODE hg38**: v46 gene annotations (feather format from AlphaGenome)
- **Cell Ontology**: OBO Foundry — CL:0000202 (IHC), CL:0000601 (OHC)

### Author Notes

This pipeline was developed to bridge the gap between population genetics (hearing loss GWAS) and mechanistic cell biology (how transcription factors wire 3D chromatin). By identifying which variants disrupt the "hardware" (sequence-encoded boundaries), we can prioritize candidates for base editing therapies in congenital deafness.

For questions or collaboration:
- Review the `haircell_workflow.py` docstrings for implementation details
- Check `media/haircell/` for example outputs
- See `help.md` for TAD analysis interpretation guide

---

## Quick Start Checklist

- [ ] Install dependencies: `make install`
- [ ] Download mouse Micro-C: `python notebooks/00_fetch_data.py`
- [ ] Set AlphaGenome API key in `.env`
- [ ] Run test: `python haircell_workflow.py ATOH1`
- [ ] Check output: `ls media/haircell/ATOH1_*`
- [ ] Interpret concordance table (Step 3)
- [ ] Identify top ISM variant (Step 5)
- [ ] Cross-reference with ClinVar / gnomAD
- [ ] Design validation experiment

**Happy 3D genome hunting! 🧬🔬👂**

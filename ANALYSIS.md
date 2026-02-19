# HiC-TAD-Library — Analysis Pipeline Documentation

This document explains every script in the repository: what it does, how the
key algorithms work, what it reads, and what it writes.

---

## Table of Contents

1. [Project Layout](#1-project-layout)
2. [Setup & Data](#2-setup--data)
3. [Core Library (`src/`)](#3-core-library-src)
   - [loaders.py](#srcloaderspy)
   - [tad_boundaries.py](#srctad_boundariespy)
   - [polymer_sim.py](#srcpolymer_simpy)
4. [Visualization Scripts](#4-visualization-scripts)
   - [visualize_tads.py](#visualize_tadspy)
   - [visualize_boundaries.py](#visualize_boundariespy)
   - [visualize_polymer_3d.py](#visualize_polymer_3dpy)
   - [visualize_alphagenome.py](#visualize_alphagenomepy)
5. [AlphaGenome Deletion Analysis](#5-alphagenome-deletion-analysis)
   - [run_mouse_deletion.py](#run_mouse_deletionpy)
   - [deletion_scan.py](#deletion_scanpy)
   - [analyze_deletion_effects.py](#analyze_deletion_effectspy)
   - [visualize_deletion_enhanced.py](#visualize_deletion_enhancedpy)
6. [Utility Scripts](#6-utility-scripts)
   - [plot_full_tad_triangles.py](#plot_full_tad_trianglespy)
   - [tal1_example_workflow.py](#tal1_example_workflowpy)
7. [Testing](#7-testing)
8. [Makefile Reference](#8-makefile-reference)

---

## 1. Project Layout

```
HiC-TAD-Library/
├── src/                        # Core reusable library
│   ├── loaders.py              # Load cooler/mcool files
│   ├── tad_boundaries.py       # TAD detection algorithms
│   ├── polymer_sim.py          # 3D polymer simulation engine
│   └── alphagenome/            # AlphaGenome API notebooks & metadata
├── tests/                      # pytest test suite (70+ tests)
│   ├── conftest.py             # Shared fixtures
│   ├── test_loaders.py
│   ├── test_tad_boundaries.py
│   └── test_polymer_sim.py
├── data/
│   ├── raw/                    # mouse_microc.mcool (not in git)
│   └── processed/
├── media/                      # All generated PNG/HTML outputs
├── externals/                  # Git submodules (alphagenome SDK, etc.)
├── notebooks/                  # Exploratory Jupyter notebooks
│── visualize_tads.py           # Hi-C heatmaps & compartment plots
├── visualize_boundaries.py     # TAD boundary analysis plots
├── visualize_polymer_3d.py     # 3D polymer simulation & rendering
├── visualize_alphagenome.py    # AlphaGenome prediction visualizations
├── run_mouse_deletion.py       # Insulator deletion comparison (AI)
├── deletion_scan.py            # Multi-position / multi-size deletion scan
├── analyze_deletion_effects.py # Before/after structural variant analysis
├── visualize_deletion_enhanced.py  # Enhanced deletion visualizations
├── plot_full_tad_triangles.py  # Full-TAD triangle plots from real Micro-C
├── tal1_example_workflow.py    # TAL1 oncogenic variant workflow example
├── ANALYSIS.md                 # This file
├── SLIDES_PLAN.md              # Slide plan for lab presentation
├── environment.yml             # Conda environment definition
├── Makefile                    # All runnable make targets
└── CLAUDE.md                   # AI assistant context file
```

---

## 2. Setup & Data

### Conda environment

```bash
conda env create -f environment.yml
conda activate hic-analysis
make install          # installs cooler/cooltools + AlphaGenome SDK
```

Python version: **3.11**. The environment includes cooler, cooltools, bioframe,
matplotlib, plotly, scipy, numpy, pandas, and pytest.

### AlphaGenome API key

Create a `.env` file in the project root:

```
ALPHA_GENOME_API_KEY=your_key_here
```

This is read by `python-dotenv` at the top of every AlphaGenome script.

### Primary data file

All visualization and boundary-detection scripts read from:

```
data/raw/mouse_microc.mcool
```

This is a multi-resolution cooler (HDF5) file from
[4DNucleome](https://data.4dnucleome.org/) containing mouse Micro-C contact
maps (mm10 genome). The default resolution used throughout is **5,000 bp (5 kb)**.
Large-region scripts (compartments) use 25 kb or 50 kb.

### Common genomic regions

| Alias | Coordinates | Context |
|-------|------------|---------|
| `Sox11_Chr12` | chr12:26,000,000–28,000,000 | Sox11 / Edward's insulator |
| `Mir9-2_Chr13` | chr13:83,500,000–84,500,000 | Mef2c / Jingyun's insulator |
| `Compartments_Chr2` | chr2:0–50,000,000 | A/B compartment example |

---

## 3. Core Library (`src/`)

### `src/loaders.py`

**What it does:** Thin wrapper that resolves the absolute path to
`data/raw/<file>` and opens a `cooler.Cooler` object at the right resolution.

```python
clr = get_cooler('mouse_microc.mcool', resolution=5000)
```

For `.mcool` files the function appends `::resolutions/<res>` to the path;
for plain `.cool` files it opens directly. Raises `FileNotFoundError` if the
file is missing so you get a clear message rather than a cryptic HDF5 error.

---

### `src/tad_boundaries.py`

The main algorithmic module. Contains five functional groups:

#### 1. Helpers

`parse_coordinates(s)` — converts `"chr12:26,000,000-28,000,000"` to
`("chr12", 26000000, 28000000)`. Strips commas so human-readable strings work.

`make_view_df(chrom, start, end)` — builds the single-row DataFrame that
cooltools functions need to know which region to operate on.

#### 2. Directionality Index (Dixon et al. 2012)

```python
di_df = compute_directionality_index(clr, "chr12:26,000,000-28,000,000")
```

For each bin **i** at resolution *r*, sums contacts in the upstream window
**A** = contacts(i − W .. i) and downstream window **B** = contacts(i .. i + W)
(default W = 500 kb). The index is:

```
E = (A + B) / 2
DI = sign(B − A) × [(A − E)² + (B − E)²] / E
```

Positive DI → more downstream contacts (inside a TAD); negative DI → more
upstream contacts. Sign changes mark TAD boundaries. Returns a DataFrame with
columns `chrom, start, end, upstream_sum, downstream_sum, DI`.

#### 3. Boundary Prominence Scoring

```python
scored = score_boundary_prominence(insulation_table, window_bp=50_000)
```

Takes an insulation score track (from cooltools) and finds local minima using
`scipy.signal.find_peaks` on the negated signal. Each minimum is assigned a
**topographic prominence** — how much the score dips below the surrounding
landscape. Boundaries are classified:

| Class | Prominence threshold |
|-------|---------------------|
| `strong` | ≥ 0.5 |
| `weak` | ≥ 0.2 |
| `sub_threshold` | < 0.2 |

The auto-threshold mode uses Li's iterative minimum cross-entropy method
(`_li_threshold`) to pick the prominence cutoff from the data itself, without
needing scikit-image.

#### 4. TAD Interval Calling

```python
tads = call_tad_intervals(boundary_df, clr, boundary_class_min='weak')
```

Pairs consecutive boundaries to define TAD intervals. Each interval records its
chromosome, start/end coordinates, length, number of bins, and the class of its
left and right boundary. Intervals longer than `max_tad_length_bp` (default
3 Mb) are discarded as likely mis-calls.

#### 5. Multi-scale Insulation & Boundary Persistence

```python
table, windows = multiscale_insulation(clr, region, n_windows=8)
persistence_df = classify_boundary_persistence(table, windows)
```

Computes the insulation score at **n** log-spaced window sizes (default 8
windows from 25 kb to 500 kb) and asks, at each genomic bin: is it a boundary
at this scale? Boundaries are then classified:

| Type | Meaning |
|------|---------|
| `constitutive` | Boundary at ≥ 75% of scales — strong, scale-invariant |
| `nested_sub` | Only at small scales — fine-grained inner structure |
| `nested_meta` | Only at large scales — coarse domain boundary |
| `mixed` | Present at some intermediate set of scales |

#### 6. Boundary Pileup

```python
pileup, n = boundary_pileup(clr, boundary_df, flank_bins=40)
```

Meta-analysis: extracts a 81×81 sub-matrix centred on each boundary and
averages them. Optionally normalises by the distance-dependent expected contact
frequency so the result shows enrichment above background. Useful for checking
that called boundaries look like real TAD boundaries (depleted diagonal,
strong flanking domains).

---

### `src/polymer_sim.py`

A self-contained 3D polymer simulation engine written in pure NumPy/SciPy
(no OpenMM, no polychrom). Converts a Hi-C contact matrix into a 3D bead-on-
a-string model.

#### Pipeline

```
Hi-C matrix
    ↓  contact_matrix_to_restraints()
Harmonic spring list [(i, j, rest_len, k), ...]
    ↓  simulate_polymer()
3D bead coordinates  (n_beads × 3)
```

#### 1. Contact matrix → restraints

Only contacts above the 70th percentile (configurable) generate springs.
Spring constant **k** scales logarithmically with contact frequency:
stronger contacts → stiffer springs. Rest length decreases with contact
strength (high frequency → shorter preferred distance), so frequently
interacting loci are pulled together.

#### 2. Backbone stiffness from insulation

`insulation_to_backbone_stiffness()` maps the 1D insulation score to per-bead
backbone spring constants. Low insulation (strong boundary) → stiff backbone
spring → a kink in the polymer chain at the TAD boundary.

#### 3. Overdamped Langevin dynamics

`simulate_polymer()` integrates:

```
Δx = (F_backbone + F_restraints + F_excluded_volume) / γ  ·  dt  +  noise
```

where γ is friction, noise is `√(2 k_B T / γ dt)` Gaussian, and
excluded-volume is a soft repulsive potential that prevents bead overlap.
No velocities are tracked (overdamped / high-friction limit, appropriate
for chromatin at genomic timescales). Runs for `n_steps` steps (default 5,000).

#### End-to-end convenience function

```python
coords = polymer_from_cooler(clr, "chr12:26,000,000-28,000,000")
# returns (n_beads, 3) array
```

---

## 4. Visualization Scripts

All visualization scripts follow the same pattern:
1. Open `data/raw/mouse_microc.mcool` at 5 kb
2. Process with `src/` modules
3. Save figures to `media/`

### `visualize_tads.py`

Produces the classic Hi-C contact map visualizations.

| Function | Output | Description |
|----------|--------|-------------|
| `plot_region()` | heatmap PNG | Log-scaled contact map for a genomic region |
| `plot_triangular_region()` | triangle PNG | 45°-rotated view, TADs appear as triangles |
| `plot_compartments()` | compartment PNG | A/B compartment eigenvector (E1) track |
| `plot_insulation()` | insulation PNG | Insulation score track with boundary calls |
| `plot_dots()` | dots PNG | Loop dot annotation overlay |
| `plot_saddle()` | saddle PNG | Compartment strength saddle plot |

Run: `python visualize_tads.py` or `make run-tads`

---

### `visualize_boundaries.py`

In-depth boundary analysis.

| Function | Output | Description |
|----------|--------|-------------|
| `plot_directionality_index()` | DI PNG | Per-bin directionality index track |
| `plot_boundary_strength()` | strength PNG | Boundary prominence ranked by strength |
| `plot_multiscale_insulation()` | multiscale PNG | Insulation heatmap across window sizes |
| `plot_tad_overlay()` | overlay PNG | Contact map with TAD intervals drawn |
| `plot_boundary_pileup()` | pileup PNG | Averaged meta-boundary contact pattern |
| `plot_combined_boundary_panel()` | combined PNG | Multi-panel summary figure |

Run: `python visualize_boundaries.py` or `make run-boundaries`

---

### `visualize_polymer_3d.py`

Converts the contact matrix to a 3D polymer model and renders it.

| Function | Output | Description |
|----------|--------|-------------|
| `run_polymer_viz()` | — | Orchestrator: runs simulation, calls renderers |
| `plot_polymer_matplotlib()` | static PNG | 3D scatter coloured by compartment/TAD |
| `plot_polymer_with_heatmap()` | side-by-side PNG | 3D model + contact map |
| `plot_dna_like_polymer_matplotlib()` | DNA-style PNG | Helical bead rendering |
| `plot_polymer_plotly()` | interactive HTML | Rotatable 3D plot in browser |
| `plot_dna_like_polymer_plotly()` | interactive HTML | DNA-style interactive plot |

Bead colouring options: A/B compartment eigenvector, insulation score, or
TAD membership. The interactive HTML files are the most useful for exploration.

Run: `python visualize_polymer_3d.py` or `make run-polymer`

---

### `visualize_alphagenome.py`

Queries the AlphaGenome API and plots the predictions.

| Function | Output | Description |
|----------|--------|-------------|
| `visualize_gene_expression()` | expression PNG | Predicted CAGE / RNA-seq tracks |
| `visualize_contact_maps()` | contact PNG | Predicted Hi-C contact map |
| `visualize_virtual_4c()` | 4C PNG | Virtual 4C from a viewpoint locus |
| `visualize_ctcf_signal()` | CTCF PNG | Predicted CTCF binding signal |

Requires `ALPHA_GENOME_API_KEY` in `.env`. Run: `python visualize_alphagenome.py`
or `make run-alphagenome`

---

## 5. AlphaGenome Deletion Analysis

These scripts use the AlphaGenome API to predict what happens to 3D chromatin
structure when a specific DNA sequence is deleted. All require the API key.

### Background: how AlphaGenome variant prediction works

AlphaGenome accepts a 1 Mb DNA sequence and predicts dozens of genomic
tracks (contact maps, CTCF binding, CAGE expression, histone marks). For
deletion analysis:

1. Predict the **wild-type (WT)** contact map for a 1 Mb window centred on the
   deletion site.
2. Simulate the deletion as a `Variant`: replace the deleted bases with a
   single `N` (effectively shortening the sequence). Predict the
   **deletion (DM)** contact map.
3. Compare WT vs DM — the difference reveals which contacts were maintained by
   the deleted sequence.

The maximum supported window size is **1,048,576 bp (2²⁰ ≈ 1 Mb)**.

---

### `run_mouse_deletion.py`

The main insulator deletion prediction script. Analyses two confirmed insulator
sites:

| Label | Coordinates | Size |
|-------|------------|------|
| Jingyun — chr13 | chr13:83,739,797–83,745,138 | ~5.3 kb |
| Edward — chr12 | chr12:27,333,532–27,336,455 | ~3 kb |

For each region × cell type, it generates **5 figure types**:

| Figure type | File pattern | Description |
|------------|--------------|-------------|
| WT \| Deletion \| Difference heatmaps | `*_<celltype>.png` | Square contact maps with gene track |
| Log₂ ratio + Virtual 4C + P(s) | `*_<celltype>_extra.png` | Quantitative measures |
| Cell-type comparison grid | `*_celltype_comparison.png` | All cell types in one figure |
| Triangle TAD view (per cell type) | `*_<celltype>_triangle.png` | Rotated 45° contact maps |
| Triangle comparison (all cell types) | `*_triangle_comparison.png` | WT/Del/Diff for each cell type |

At the end it writes `analysis_report.html` — a self-contained HTML report
with all figures.

**Cell types used:**

| ID | Name |
|----|------|
| `CL:0000207` | Olfactory receptor cell |
| `EFO:0004038` | Mouse embryonic stem cell (PI recommendation) |

Run: `python run_mouse_deletion.py` or `make run-deletion`

---

### `deletion_scan.py`

Tests **12 evenly-spaced deletion positions** across the 1 Mb window to ask:
*is the actual insulator site special, or would any deletion produce a similar
effect?*

#### Usage

```bash
python deletion_scan.py edward                         # original ~3 kb size
python deletion_scan.py jingyun                        # original ~5.3 kb size
python deletion_scan.py edward  --del-sizes 10 40 80   # 10/40/80 kb
python deletion_scan.py jingyun --del-sizes 10 40 80
```

#### How it works

1. **Predict WT once** — reused across all deletion sizes.
2. **Build scan positions** — `np.linspace` across the usable window, then
   snap the nearest position to the actual insulator centre.
3. **For each position × size** — call `predict_variant` with a deletion of
   that size centred at that position.
4. **Compute three metrics** for each (WT, DM) pair:

| Metric | Formula | Interpretation |
|--------|---------|----------------|
| `global_abs` | mean \|DM − WT\| | Total contact reorganisation |
| `cross_gain` | mean(DM[left, right]) − mean(WT[left, right]) | Cross-boundary contact change; positive = TADs merging |
| `ins_weakening` | mean(DM[diamond]) − mean(WT[diamond]) | Change in contacts in a 50 kb diamond centred on the deletion; positive = boundary weakened |

5. **Rank positions** by `global_abs` and generate four figures:

| Figure | File | Description |
|--------|------|-------------|
| Sensitivity profile | `deletion_scan_<region>_<size>kb_sensitivity.png` | Bar chart of all 3 metrics |
| Ranked summary | `deletion_scan_<region>_<size>kb_ranked_summary.png` | WT triangle + impact ranking bars |
| Triangle gallery | `deletion_scan_<region>_<size>kb_triangle_gallery.png` | Difference maps for selected sites |
| Diff montage | `deletion_scan_<region>_<size>kb_montage.png` | All 12 sites on a shared scale |

When multiple `--del-sizes` are given, a fifth figure is also produced:

| Figure | File | Description |
|--------|------|-------------|
| Cross-size comparison | `deletion_scan_<region>_cross_size_comparison.png` | Heatmaps (position × size) + rank trend |

Finally, a self-contained HTML report is written:
`deletion_scan_<region>_<sizes>_report.html`

#### Key results

| Region | Original size (~3–5 kb) rank | 10 kb rank | 40 kb rank | 80 kb rank |
|--------|:---:|:---:|:---:|:---:|
| Edward chr12 | #7/12 | **#1/12** | **#1/12** | #4/12 |
| Jingyun chr13 | #7/12 | #7/12 | #7/12 | #7/12 |

Edward's chr12 insulator becomes rank #1 with 10–40 kb deletions (insulation
weakening +0.46, ~10× any other site), confirming it is a genuine functional
boundary. Jingyun's chr13 site stays at #7 at all sizes, suggesting a
different or more redundant architectural role in Mouse ESC.

---

### `analyze_deletion_effects.py`

Computes structural variant effects on the contact map for a configurable
set of deletions. Simpler than `deletion_scan.py` — focuses on visualising
before/after for a single deletion rather than scanning many positions.

`visualize_deletion_effect()` — produces a three-panel figure: WT, DM,
and difference map (RdBu colour scale).

Run: `python analyze_deletion_effects.py` or `make run-deletion-effects`

---

### `visualize_deletion_enhanced.py`

Extended version of the deletion visualisation with additional tracks and
annotation overlays (gene models, CTCF peaks, insulation score). Intended
for producing publication-quality panels.

Run: `python visualize_deletion_enhanced.py` or `make run-deletion-enhanced`

---

## 6. Utility Scripts

### `plot_full_tad_triangles.py`

Plots Edward's and Jingyun's complete 2-TAD pair regions as triangle contact
maps using **real Mouse Micro-C data** (not AlphaGenome predictions).

| Region | Coordinates | Bins |
|--------|------------|------|
| Edward chr12 | chr12:26,440,002–28,560,000 | 424 × 424 |
| Jingyun chr13 | chr13:81,760,002–85,200,000 | 688 × 688 |

The triangle rotation formula maps matrix element `[i, j]` (j ≥ i) to:

```
x = i + j    (genomic midpoint of the contact pair)
y = j − i    (genomic distance / 2)
```

producing a (n, 2n) array displayed with `origin='lower'` so the diagonal
sits at the bottom and TADs appear as dark upward-pointing triangles.

Outputs:
- `media/Edward_chr12_full_tad_triangle.png`
- `media/Jingyun_chr13_full_tad_triangle.png`
- `media/full_tad_triangles_combined.png`

Run: `python plot_full_tad_triangles.py` or `make run-tad-triangles`

---

### `tal1_example_workflow.py`

Demonstrates the AlphaGenome variant effect prediction workflow using the
TAL1 oncogenic translocation (chr1:47,595,955 T→A) as an example.

Steps:
1. Load mm10 GENCODE M23 gene annotations with `load_gtf_and_extractor()`
2. Define the TAL1 variant(s) with `oncogenic_tal1_variants()`
3. Query AlphaGenome for predicted RNA expression, accessibility, and
   histone modifications at each variant position
4. Generate multi-track comparison figures

This is a worked example / tutorial for how to extend the deletion analysis
to point mutations and other structural variants.

Run: `python tal1_example_workflow.py` or `make run-tal1`

---

## 7. Testing

```bash
make test              # run all tests
make test-unit         # fast unit tests only (no data needed)
make test-coverage     # generate HTML coverage report in htmlcov/
make test-verbose      # detailed per-test output
```

### Structure

```
tests/
├── conftest.py              # Shared fixtures (mock matrices, coolers, etc.)
├── test_loaders.py          # 10 tests  — src/loaders.py
├── test_tad_boundaries.py   # 35+ tests — src/tad_boundaries.py
└── test_polymer_sim.py      # 25+ tests — src/polymer_sim.py
```

### Test markers

| Marker | Meaning |
|--------|---------|
| `@pytest.mark.unit` | Fast, no external files needed |
| `@pytest.mark.integration` | May use mock coolers |
| `@pytest.mark.slow` | Longer simulations |
| `@pytest.mark.requires_data` | Needs real mcool file |

All unit tests use mock data from `conftest.py` — a synthetic 100×100
contact matrix with two TAD blocks and a boundary gap, plus corresponding
mock insulation scores and mock cooler objects.

---

## 8. Makefile Reference

See the Makefile for the full target list. Below is a quick reference:

```bash
# ── Setup ──────────────────────────────────────────────────────────
make install                   # install all dependencies
make setup                     # install + run core visualizations

# ── Core visualizations (use real Micro-C data) ────────────────────
make run-tads                  # visualize_tads.py
make run-boundaries            # visualize_boundaries.py
make run-polymer               # visualize_polymer_3d.py
make run-alphagenome           # visualize_alphagenome.py
make run                       # tads + boundaries + alphagenome

# ── AlphaGenome deletion analysis (require API key) ────────────────
make run-deletion              # run_mouse_deletion.py (both regions)
make run-deletion-effects      # analyze_deletion_effects.py
make run-deletion-enhanced     # visualize_deletion_enhanced.py

# ── Deletion sensitivity scans (original insulator size) ───────────
make scan-edward               # edward, ~3 kb deletions
make scan-jingyun              # jingyun, ~5.3 kb deletions
make scan-all                  # both regions, original size

# ── Deletion sensitivity scans (10 / 40 / 80 kb) ──────────────────
make scan-large-edward         # edward, 10 40 80 kb
make scan-large-jingyun        # jingyun, 10 40 80 kb
make scan-large                # both regions, large sizes

# ── Utility scripts ────────────────────────────────────────────────
make run-tad-triangles         # real Micro-C full-TAD triangles
make run-tal1                  # TAL1 variant workflow example

# ── Run everything ─────────────────────────────────────────────────
make run-all                   # all analysis (no deletion scans)
make analysis                  # core viz + both deletion scans

# ── Testing ────────────────────────────────────────────────────────
make test                      # all tests
make test-unit                 # unit tests only
make test-coverage             # with HTML coverage report
make test-verbose              # detailed output

# ── Housekeeping ───────────────────────────────────────────────────
make clean                     # remove media/*.png and test artefacts
make clean-html                # remove generated HTML reports
make clean-all                 # media/ + HTML + htmlcov/ + caches
make help                      # print this reference
```

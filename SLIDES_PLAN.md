# Presentation Slide Plan
## "Predicting Insulator Deletion Effects with AlphaGenome"
### Informal lab presentation — 1–2 slides

---

## Slide 1: What I'm Doing

**Title:** Using AI to Predict What Happens When You Delete a TAD Insulator

**Talking points:**
- 3D genome is organized into blocks called **TADs** (Topologically Associating Domains)
- TADs are separated by **insulators** — short regions (~5 kb) that act like walls
- **Prediction:** deleting an insulator should merge the two neighboring TADs
- **AlphaGenome** (Google DeepMind) is an AI model that predicts chromatin contact maps from raw DNA sequence
- Feed it the DNA sequence → get a predicted Hi-C-style contact map
- Feed it the same sequence *with the deletion* → see if the TADs merge

**Visual ideas (from `media/` folder):**

*Triangle TAD view (most intuitive for an audience):*
- `mouse_deletion_chr13_..._triangle_comparison.png` — WT | Deletion | Difference
  as rotated-45° triangles. TADs appear as dark upward triangles; the insulator
  shows as a valley between them. Cyan lines mark the deletion site.
- `mouse_deletion_chr12_..._triangle_comparison.png` — same for Ed's region

*Square heatmap (classic Hi-C view):*
- `mouse_deletion_chr13_..._celltype_comparison.png`
- `mouse_deletion_chr12_..._celltype_comparison.png`

*Full interactive report:* `analysis_report.html` (open in browser)

---

## Slide 2: The Two Experiments

**Title:** Two Insulator Deletions — Jingyun's & Ed's Regions (mm10)

**Jingyun's region (Chr13):**
- Insulator: `chr13:83,739,797 – 83,745,138` (~5.3 kb deletion)
- Full 2-TAD window: `chr13:81,760,002 – 85,200,000` (3.4 Mb span)
- Prediction: merging of the two TADs flanking this insulator

**Ed's region (Chr12):**
- Insulator: `chr12:27,333,532 – 27,336,455` (~3 kb deletion)
- Full 2-TAD window: `chr12:26,440,002 – 28,560,000` (2.1 Mb span)
- Prediction: merging of the two TADs flanking this insulator

**Cell types used in the model:**
- `CL:0000207` — Olfactory receptor cell (cranial placode origin, closest proxy)
- `EFO:0004038` — Mouse embryonic stem cell (suggested by PI)
- Note: No inner ear / cochlear / otic placode cell type exists in AlphaGenome's
  training data. Only 8 mouse contact map cell types are available.
  `CL:0000589` (cochlear inner hair cell) is available for CAGE tracks but not contact maps.

**Note on window size:** AlphaGenome supports a maximum of 1,048,576 bp (1 Mb).
The full 2-TAD windows (3.4 Mb / 2.1 Mb) are too large for the model API, so
predictions are made over a 1 Mb window centred on the insulator deletion site.

**Visual ideas:**
- Triangle comparison plot for each region (most slide-friendly)
- Virtual 4C + P(s) curve from `_extra.png` files
- Full HTML report: `analysis_report.html` (open in browser for all figures)

---

## Code & GitHub

The analysis lives in [`run_mouse_deletion.py`](run_mouse_deletion.py):
- Loads mm10 GENCODE M23 gene annotations
- Queries AlphaGenome API for each cell type × each region
- Generates 5 figure types per run:
  1. WT | Deletion | Difference heatmaps (square, with gene track)
  2. Log₂ ratio map + Virtual 4C + P(s) curve
  3. Cell-type comparison grid (square heatmaps)
  4. **Triangle TAD view** — per cell type (rotated 45°, TADs as triangles)
  5. **Triangle comparison** — all cell types × WT | Deletion | Difference

All outputs saved to `media/` as PNG files.
Full HTML report with all figures: [`analysis_report.html`](analysis_report.html)

### How the triangle plot works
```
_rotate45(matrix):
  For each upper-triangle element matrix[i,j] (j >= i):
    x_rotated = i + j      → genomic midpoint (horizontal axis)
    y_rotated = j - i      → genomic distance / 2 (vertical axis)
  Result: (n, 2n) array displayed with origin='lower'
  → diagonal at bottom, TADs as dark upward-pointing triangles
```

---

## Ontology Notes (for reference)

| Term | Name | In AlphaGenome contact maps? |
|------|------|------------------------------|
| `CL:0000207` | Olfactory receptor cell | Yes |
| `EFO:0004038` | Mouse embryonic stem cell | Yes (suggested by PI) |
| `CL:0000589` | Cochlear inner hair cell | CAGE only, not contact maps |
| `UBERON:0001846` | Internal ear | CAGE only, not contact maps |
| `UBERON:0003069` | Otic placode | Not in training data |
| `GO:1905040` | Otic placode development | Not in training data |

---

## Key Result to Show

> AlphaGenome predicts that deleting a ~5 kb insulator sequence causes the two
> flanking TADs to partially merge — visible as new off-diagonal contacts in the
> heatmap and a gain of interactions at long genomic distances in the P(s) curve.

This is a *computational prediction* to complement the wet-lab deletion experiments.

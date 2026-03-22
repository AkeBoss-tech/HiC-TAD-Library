# Synthesis Experiments — In-silico Chromatin Engineering

Three quick validation experiments demonstrating AI-driven chromatin boundary engineering
at the **Unc5b CTCF cluster** (mm10 chr10:60,240,000–61,288,576).

## Prerequisites

- Conda environment: `hic-analysis` (Python 3.11)
- `ALPHA_GENOME_API_KEY` in `.env`
- Internet access (AlphaGenome API calls)

## How to Run

```bash
conda activate hic-analysis

# Run all three experiments + generate report
make run-synthesis

# Or run individually
make run-synthesis-exp1     # ~15 s  — 1 API call (WT)
make run-synthesis-exp2     # ~45 s  — 3 API calls
make run-synthesis-exp3     # ~75 s  — 5 API calls
make synthesis-report       # instant — generates synthesis_report.html
```

## Outputs

| File | Description |
|------|-------------|
| `media/synthesis_exp1_baseline_unc5b.png` | WT contact map · CTCF signal · insulation |
| `media/synthesis_exp2_contact_panels.png` | 2×4 contact panels + difference maps |
| `media/synthesis_exp2_delta_insulation.png` | Δ insulation bar chart |
| `media/synthesis_exp3_synthetic_insulator.png` | Insulation curves + best design comparison |
| `synthesis_report.html` | Full HTML report |
| `data/processed/synthesis_unc5b_wt_cache.npz` | WT cache (avoids re-calling API) |

## Experiment Summary

| # | Script | API calls | What it tests |
|---|--------|-----------|---------------|
| 1 | `exp1_baseline.py` | 1 | WT CTCF cluster prediction |
| 2 | `exp2_ctcf_deletions.py` | 3 | Insulation loss per CTCF site removed |
| 3 | `exp3_synthetic_insulator.py` | 5 | Convergent CTCF pair insertions |

**Total:** 9 API calls, ~2 min end-to-end.

## Design Details

**Locus:** Unc5b CTCF cluster — 4–6 convergent CTCF sites validated in the 2025 AlphaGenome paper.
**Cell type:** `EFO:0004038` (Mouse ESC) — consistent with all existing deletion experiments.
**CTCF motifs:** `CCGCGAGGCGCAGG` (FWD) / `CCTGCGCCTCGCGG` (REV, convergent pair).
**WT cache:** After Exp 1, subsequent experiments load from `data/processed/synthesis_unc5b_wt_cache.npz`
to avoid redundant API calls.

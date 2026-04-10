# HG-DT v1.1 roadmap

This document tracks the upgrade path from the current **HG-DT** Streamlit prototype (`app.py` + `src/hg_dt/`) toward the **v1.1 specification** (local reference data, flexible edits, full DNA→RNA→protein, richer causality and visualization).

## Phase 0 — Setup and safety net

- [x] Branch `v1.1-full-spec` for integration work.
- [x] Dependency notes: `requirements.txt` (pip) complements `environment.yml` (conda).
- [x] `pip install -r requirements.txt` and HG-DT unit tests.
- [ ] Run `make hg-dt-ui` after each milestone to confirm the wizard still runs.

## Phase 1 — De-hardcode and polish (short term)

- [x] User-controlled **insertion** sequence (validated length and ACGT); **SNP** alternate base selector.
- [x] Centralized sequence logic in `src/hg_dt/data/sequence_fetcher.py` (`fetch_ref_mut_sequences`, `build_alpha_window`, `compute_mut_shift`).
- [x] Wire **`predict_isoforms`** + **`compare_translation`** into the protein path (CDS fallback retained).
- [x] Trajectory caption aligned with **H1 + GM12878** contact maps.
- [x] ESMFold **~400 AA** limit surfaced in the UI with truncation warning.

## Phase 2 — Local reference and context builder (≈1 week)

- [x] `make download-references` — `scripts/download_hgdt_references.sh` (hg38, GENCODE 49 GTF, SCREEN ELS BED). See `docs/hg-dt/REFERENCES.md`.
- [x] Local **hg38 FASTA** drives AlphaGenome sequence windows when `data/references/hg38.fa` is present (`fetch_ref_mut_sequences(..., fasta_path=...)` + `reference_paths.py`).
- [ ] Wire `ReferenceContextBuilder` + GENCODE into Step 2 / isoform path (optional UCSC fallback remains for zero-setup demos).

## Phase 3 — DNA → RNA → protein completeness (≈1 week)

- [ ] Tie isoform-level interpretation to AlphaGenome expression tracks where available.
- [ ] Optional second structure backend (e.g. AlphaFold 3 / BioNeMo) alongside ESMFold for long proteins.

## Phase 4 — Visualization and interpretation (≈3–5 days)

- [ ] Streamlit tabs / `render_full_report`-style layout (`src/hg_dt/viz/dashboard.py` candidate).
- [ ] pyGenomeTracks-based 1D browser; optional NGLview/py3Dmol for proteins; trajectory GIF export.

## Phase 5 — Research-grade polish (≈1–2 weeks)

- [ ] Uncertainty on deltas; stronger attribution using full cCRE + gene overlap.
- [ ] Broader cell-type discovery from AlphaGenome ontology; OOD / locus quality flags.
- [ ] CLI batch runner (e.g. `hg-dt run --locus …`).

---

*Last updated: April 2026.*

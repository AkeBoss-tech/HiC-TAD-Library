# Makefile for HiC-TAD-Library
# Run `make help` to see all available targets.
#
# Prerequisites:
#   conda activate hic-analysis   (or ensure the right python is on PATH)
#   ALPHA_GENOME_API_KEY=...      set in .env for AlphaGenome targets

PYTHON  ?= python
PYTEST  ?= pytest
PIP     ?= pip

# Deletion sizes used for the large-deletion scans
LARGE_SIZES ?= 10 40 80

.PHONY: help \
        install install-alphagenome \
        run run-tads run-boundaries run-polymer run-alphagenome \
        run-deletion run-deletion-effects run-deletion-enhanced \
        scan-edward scan-jingyun scan-all \
        scan-large-edward scan-large-jingyun scan-large \
        run-tad-triangles run-tal1 \
        run-compartments run-loops run-sv-cnv run-variants \
        run-hichip run-capture-hic run-phasing run-dovetail \
        run-all analysis \
        test test-unit test-integration test-coverage test-verbose \
        clean clean-html clean-all \
        setup docs

# ── Help ──────────────────────────────────────────────────────────────────────
help:
	@echo ""
	@echo "HiC-TAD-Library — available make targets"
	@echo "========================================="
	@echo ""
	@echo "Setup"
	@echo "  make install              Install all Python dependencies"
	@echo "  make setup                Install + run core visualizations"
	@echo ""
	@echo "Core visualizations  (use real Mouse Micro-C data)"
	@echo "  make run-tads             visualize_tads.py"
	@echo "  make run-boundaries       visualize_boundaries.py"
	@echo "  make run-polymer          visualize_polymer_3d.py"
	@echo "  make run-alphagenome      visualize_alphagenome.py"
	@echo "  make run                  tads + boundaries + alphagenome"
	@echo ""
	@echo "AlphaGenome deletion analysis  (require API key in .env)"
	@echo "  make run-deletion         run_mouse_deletion.py  (both regions)"
	@echo "  make run-deletion-effects analyze_deletion_effects.py"
	@echo "  make run-deletion-enhanced visualize_deletion_enhanced.py"
	@echo ""
	@echo "Deletion sensitivity scans  (original insulator size)"
	@echo "  make scan-edward          Edward chr12, ~3 kb deletions"
	@echo "  make scan-jingyun         Jingyun chr13, ~5.3 kb deletions"
	@echo "  make scan-all             Both regions, original size"
	@echo ""
	@echo "Deletion sensitivity scans  (10 / 40 / 80 kb)"
	@echo "  make scan-large-edward    Edward chr12, $(LARGE_SIZES) kb"
	@echo "  make scan-large-jingyun   Jingyun chr13, $(LARGE_SIZES) kb"
	@echo "  make scan-large           Both regions, large sizes"
	@echo ""
	@echo "Utility scripts"
	@echo "  make run-tad-triangles    Real Micro-C full-TAD triangle plots"
	@echo "  make run-tal1             TAL1 oncogenic variant workflow"
	@echo "  make docs                 Generate HTML documentation from Markdown files"
	@echo ""
	@echo "Dovetail analysis suite"
	@echo "  make run-compartments     A/B compartment analysis (real mcool)"
	@echo "  make run-loops            Chromatin loop calling (real mcool)"
	@echo "  make run-sv-cnv           SV and CNV detection (real mcool + optional BAM)"
	@echo "  make run-variants         SNV/indel calling (BAM; demo mode if absent)"
	@echo "  make run-hichip           HiChIP loop analysis (real mcool; demo peaks if absent)"
	@echo "  make run-capture-hic      Capture Hi-C interaction calling (demo baits if absent)"
	@echo "  make run-phasing          Hi-C haplotype phasing (BAM+VCF; demo if absent)"
	@echo "  make run-dovetail         All seven Dovetail analyses"
	@echo ""
	@echo "Combined targets"
	@echo "  make run-all              All visualizations (no deletion scans)"
	@echo "  make analysis             Core viz + both original-size scans"
	@echo ""
	@echo "Testing"
	@echo "  make test                 Run all tests"
	@echo "  make test-unit            Unit tests only  (no data required)"
	@echo "  make test-integration     Integration tests"
	@echo "  make test-coverage        Tests + HTML coverage report"
	@echo "  make test-verbose         Detailed per-test output"
	@echo ""
	@echo "Housekeeping"
	@echo "  make clean                Remove media/*.png and test artefacts"
	@echo "  make clean-html           Remove generated HTML reports"
	@echo "  make clean-all            media/ + HTML + htmlcov/ + caches"
	@echo ""

# ── Installation ──────────────────────────────────────────────────────────────
install-alphagenome:
	$(PIP) install -e ./externals/alphagenome
	$(PIP) install python-dotenv

install: install-alphagenome
	$(PIP) install cooler cooltools bioframe matplotlib pandas numpy scipy plotly

# ── Core visualizations ───────────────────────────────────────────────────────
run-tads:
	$(PYTHON) visualize_tads.py

run-boundaries:
	$(PYTHON) visualize_boundaries.py

run-polymer:
	$(PYTHON) visualize_polymer_3d.py

run-alphagenome:
	$(PYTHON) visualize_alphagenome.py

# Original `make run` target — kept for backward compatibility
run: run-tads run-boundaries run-alphagenome

# ── AlphaGenome deletion analysis ─────────────────────────────────────────────
run-deletion:
	$(PYTHON) run_mouse_deletion.py

run-deletion-effects:
	$(PYTHON) analyze_deletion_effects.py

run-deletion-enhanced:
	$(PYTHON) visualize_deletion_enhanced.py

# ── Deletion sensitivity scans (original insulator size) ──────────────────────
scan-edward:
	$(PYTHON) deletion_scan.py edward

scan-jingyun:
	$(PYTHON) deletion_scan.py jingyun

scan-all: scan-edward scan-jingyun

# ── Deletion sensitivity scans (10 / 40 / 80 kb) ─────────────────────────────
scan-large-edward:
	$(PYTHON) deletion_scan.py edward --del-sizes $(LARGE_SIZES)

scan-large-jingyun:
	$(PYTHON) deletion_scan.py jingyun --del-sizes $(LARGE_SIZES)

scan-large: scan-large-edward scan-large-jingyun

# ── Utility scripts ───────────────────────────────────────────────────────────
run-tad-triangles:
	$(PYTHON) plot_full_tad_triangles.py

run-tal1:
	$(PYTHON) tal1_example_workflow.py

docs:
	$(PYTHON) scripts/build_docs.py

# ── Dovetail analysis suite ───────────────────────────────────────────────────
run-compartments:
	$(PYTHON) visualize_compartments.py

run-loops:
	$(PYTHON) visualize_loops.py

run-sv-cnv:
	$(PYTHON) analyze_sv_cnv.py

run-variants:
	$(PYTHON) analyze_variants.py

run-hichip:
	$(PYTHON) analyze_hichip.py

run-capture-hic:
	$(PYTHON) analyze_capture_hic.py

run-phasing:
	$(PYTHON) analyze_phasing.py

# Run all Dovetail analyses (compartments + loops use real data;
# variants/hichip/capture/phasing use demo data if BAM/VCF/BED absent)
run-dovetail: run-compartments run-loops run-sv-cnv run-variants \
              run-hichip run-capture-hic run-phasing

# ── Combined targets ──────────────────────────────────────────────────────────

# All visualizations that use real Micro-C data (no API key needed)
run-all: run-tads run-boundaries run-polymer run-alphagenome run-tad-triangles

# Full analysis pipeline: core visualizations + deletion scans
analysis: run-tads run-boundaries run-deletion scan-all

# ── Testing ───────────────────────────────────────────────────────────────────
test:
	$(PYTEST)

test-unit:
	$(PYTEST) -m unit

test-integration:
	$(PYTEST) -m integration

test-coverage:
	$(PYTEST) --cov=src --cov-report=html --cov-report=term

test-verbose:
	$(PYTEST) -vv

# ── Housekeeping ──────────────────────────────────────────────────────────────
clean:
	rm -f media/*.png
	rm -rf htmlcov/ .coverage
	find . -type d -name __pycache__ -exec rm -rf {} + 2>/dev/null || true
	find . -type d -name .pytest_cache -exec rm -rf {} + 2>/dev/null || true

clean-html:
	rm -f analysis_report.html
	rm -f deletion_scan_*.html
	rm -f README.html ANALYSIS.html TESTING.html notebooks.html

clean-all: clean clean-html
	rm -f media/*.html

# ── Setup (install + first run) ───────────────────────────────────────────────
setup: install run

# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

HiC-TAD-Library is a Python library for visualizing and analyzing Hi-C chromatin contact data, focusing on Topologically Associating Domains (TADs), compartments, and insulation scores. The project integrates AlphaGenome AI for high-resolution genomic predictions and features a custom 3D polymer simulation engine.

## Build & Development Commands

### Setup
```bash
make setup              # Full setup: install dependencies + run all visualizations
make install            # Install all dependencies (includes AlphaGenome)
make install-alphagenome  # Install AlphaGenome package from externals/
```

### Running Visualizations
```bash
make run                # Run all visualization scripts (tads, boundaries, alphagenome)
python visualize_tads.py           # TAD and compartment visualizations
python visualize_boundaries.py     # Boundary analysis (strength, pileup, combined)
python visualize_polymer_3d.py     # 3D polymer simulations from contact matrices
python visualize_alphagenome.py    # AlphaGenome predictions
python tal1_example_workflow.py    # TAL1 variant analysis workflow
```

### Data Management
```bash
make clean              # Remove all generated PNG files from media/
```

### Environment
- Conda environment: `hic-analysis` (defined in environment.yml)
- Python version: 3.11
- AlphaGenome API key required in `.env` as `ALPHA_GENOME_API_KEY`

## Architecture

### Core Modules (`src/`)

**`src/tad_boundaries.py`** - Comprehensive TAD boundary detection toolkit:
- `compute_directionality_index()` - Dixon et al. directionality index for TAD calling
- `score_boundary_prominence()` - Ranks boundaries by insulation depth
- `call_tad_intervals()` - Converts boundary positions to TAD intervals
- `compute_pileup_around_boundaries()` - Meta-analysis of boundary contact patterns
- Supports multi-scale insulation analysis (multiple window sizes)

**`src/polymer_sim.py`** - Custom 3D polymer simulation engine:
- `contact_matrix_to_restraints()` - Converts Hi-C matrices to harmonic spring restraints
- `simulate_polymer()` - Overdamped Langevin dynamics simulation (numpy/scipy only, no external engines)
- No dependency on polychrom or OpenMM - pure Python implementation
- Generates 3D bead coordinates from contact frequencies for visualization

**`src/alphagenome/`** - AlphaGenome SDK integration:
- AlphaGenome API wrapper for genomic predictions
- Variant effect prediction (RNA, accessibility, histone modifications)
- Predicted contact maps, CTCF binding, Virtual 4C

**`src/loaders.py`** - Data loading utilities for cooler/mcool files

### Visualization Scripts (root level)

All visualization scripts follow the pattern:
1. Load Hi-C data from `data/raw/mouse_microc.mcool` at 5kb resolution
2. Process using `src/` modules
3. Generate figures in `media/` directory

**Common regions** (mm10 mouse genome):
- `Sox11_Chr12`: "chr12:26,000,000-28,000,000"
- `Mir9-2_Chr13`: "chr13:83,500,000-84,500,000"
- `Compartments_Chr2`: "chr2:0-50,000,000"

### Data Pipeline

**Data sources**: [4DNucleome](https://data.4dnucleome.org/)
- Primary dataset: Mouse Micro-C (`data/raw/mouse_microc.mcool`)
- Resolution: 5000 bp (5kb) for TAD analysis, configurable in scripts
- Format: Multi-resolution cooler (.mcool) files

**Data directories**:
- `data/raw/` - Raw Hi-C/Micro-C data files (.mcool, .cool)
- `data/processed/` - Processed analysis outputs
- `data/external/` - External reference datasets

**Output**: All generated visualizations saved to `media/` as PNG files

## Testing

### Running Tests

The project includes a comprehensive test suite with 70+ tests:

```bash
make test              # Run all tests
make test-coverage     # Run with coverage report
make test-unit         # Run only unit tests (fast, no data required)
make test-verbose      # Run with detailed output
```

### Test Organization

- **`tests/test_loaders.py`** - Data loading utilities (10 tests)
- **`tests/test_tad_boundaries.py`** - TAD boundary detection (35+ tests)
- **`tests/test_polymer_sim.py`** - 3D polymer simulation (25+ tests)
- **`tests/conftest.py`** - Shared fixtures and mock data

### Test Markers

- `@pytest.mark.unit` - Fast unit tests (no external data)
- `@pytest.mark.integration` - Integration tests (may need data files)
- `@pytest.mark.slow` - Longer-running tests
- `@pytest.mark.requires_data` - Requires actual Hi-C data files

### Coverage Reports

```bash
make test-coverage        # Generate HTML + terminal coverage report
open htmlcov/index.html   # View detailed coverage report
```

**See [TESTING.md](TESTING.md) for comprehensive testing documentation.**

### External Dependencies (`externals/`)

Git submodules for related analysis tools:
- `alphagenome/` - AlphaGenome SDK (pip install -e)
- `GAM-analysis/` - Genome Architecture Mapping benchmarking
- `loop-analysis/` - Chromatin loop analysis methods
- `methods-benchmarking/` - TAD calling method comparisons

## Key Concepts

**TAD Analysis Pipeline**:
1. Compute insulation scores at multiple window sizes (25kb, 50kb, 100kb typical)
2. Detect boundaries as local minima in insulation track
3. Score boundary prominence (strength) based on insulation depth
4. Call TAD intervals from boundary positions
5. Classify boundaries as strong/weak based on prominence thresholds

**3D Polymer Simulation Workflow**:
1. Extract balanced contact matrix from cooler file
2. Convert high-frequency contacts to harmonic restraints (threshold ~70th percentile)
3. Run Langevin dynamics with connectivity constraints
4. Color beads by: A/B compartment (E1 eigenvector), insulation score, or TAD membership
5. Render with matplotlib (static) or plotly (interactive HTML)

**Coordinate Parsing**: Use `parse_coordinates()` to handle strings like "chr12:26,000,000-28,000,000"

## Important Notes

- **File paths**: Mac-specific absolute paths hardcoded in visualization scripts (update `FILE_PATH` variable)
- **Memory**: 5kb resolution is high-resolution; use 10kb for larger regions or memory constraints
- **AlphaGenome**: Requires valid API key and internet connectivity
- **Notebooks**: Exploratory analysis in `notebooks/` (data fetching, TAD exploration, compartments)
- **Git submodules**: Run `git submodule update --init --recursive` to populate `externals/`

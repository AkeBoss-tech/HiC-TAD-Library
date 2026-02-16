# Testing Guide for HiC-TAD-Library

## Overview

The HiC-TAD-Library includes a comprehensive test suite with **70+ tests** covering all core functionality:

- **test_loaders.py** (10 tests) - Data loading utilities
- **test_tad_boundaries.py** (35+ tests) - TAD boundary detection and analysis
- **test_polymer_sim.py** (25+ tests) - 3D polymer simulation

## Quick Start

```bash
# Run all tests
make test

# Run with coverage
make test-coverage

# Run only unit tests (fast)
make test-unit

# Run verbose output
make test-verbose
```

## Installation

Install test dependencies:

```bash
# Update conda environment with test packages
conda env update -f environment.yml

# Or install manually
conda install pytest pytest-cov pytest-mock
```

## Test Commands Reference

### Makefile Commands

| Command | Description |
|---------|-------------|
| `make test` | Run all tests |
| `make test-unit` | Run only unit tests (fast, no data required) |
| `make test-integration` | Run integration tests (may require data files) |
| `make test-coverage` | Run tests with coverage report (HTML + terminal) |
| `make test-verbose` | Run tests with detailed output |

### Direct Pytest Commands

```bash
# Run all tests
pytest

# Run specific test file
pytest tests/test_loaders.py
pytest tests/test_tad_boundaries.py
pytest tests/test_polymer_sim.py

# Run specific test class
pytest tests/test_tad_boundaries.py::TestBoundaryProminence

# Run specific test function
pytest tests/test_loaders.py::TestGetCooler::test_get_cooler_default_resolution

# Run with pattern matching
pytest -k "boundary"  # All tests with "boundary" in name
pytest -k "not slow"  # Exclude slow tests

# Run with markers
pytest -m unit         # Only unit tests
pytest -m integration  # Only integration tests
pytest -m "not slow"   # Exclude slow tests

# Verbose output options
pytest -v              # Verbose
pytest -vv             # Very verbose
pytest -s              # Show print statements

# Debug options
pytest --tb=short      # Short traceback
pytest --tb=long       # Long traceback
pytest --pdb           # Drop into debugger on failure
pytest --lf            # Run last failed tests
pytest --ff            # Run failures first, then rest

# Coverage options (requires pytest-cov)
pytest --cov=src                           # Basic coverage
pytest --cov=src --cov-report=html         # HTML report
pytest --cov=src --cov-report=term-missing # Show missing lines
pytest --cov=src --cov-branch              # Branch coverage
```

## Test Structure

### Directory Layout

```
tests/
├── __init__.py           # Package marker
├── README.md             # Detailed test documentation
├── conftest.py           # Shared fixtures and test data
├── test_loaders.py       # Tests for src/loaders.py
├── test_tad_boundaries.py # Tests for src/tad_boundaries.py
└── test_polymer_sim.py   # Tests for src/polymer_sim.py
```

### Test Markers

Tests are organized with pytest markers for selective execution:

- **`@pytest.mark.unit`** - Fast unit tests, no external dependencies
- **`@pytest.mark.integration`** - Integration tests, may need data files
- **`@pytest.mark.slow`** - Tests that take longer to run
- **`@pytest.mark.requires_data`** - Tests requiring actual Hi-C data files

## Coverage Reports

### Generate Coverage Report

```bash
# Terminal + HTML report
make test-coverage

# Or with pytest directly
pytest --cov=src --cov-report=html --cov-report=term
```

### View HTML Report

```bash
# Open in browser
open htmlcov/index.html      # macOS
xdg-open htmlcov/index.html  # Linux
start htmlcov/index.html     # Windows
```

The HTML report shows:
- Line-by-line coverage for each module
- Highlighted covered/uncovered lines
- Branch coverage statistics
- Summary by module

## Test Modules

### test_loaders.py

Tests data loading functionality:

✅ Path construction for .mcool and .cool files
✅ Resolution URI formatting
✅ Error handling for missing files
✅ Default parameter values
✅ Integration with real data files (when available)

**Example:**
```bash
pytest tests/test_loaders.py -v
```

### test_tad_boundaries.py

Comprehensive tests for TAD analysis:

✅ Coordinate parsing (`parse_coordinates`)
✅ Directionality Index computation
✅ Boundary prominence scoring
✅ TAD interval calling
✅ Multi-scale insulation analysis
✅ Boundary persistence classification
✅ Boundary pileup analysis
✅ Edge cases and error handling

**Example:**
```bash
pytest tests/test_tad_boundaries.py::TestBoundaryProminence -v
```

### test_polymer_sim.py

Tests 3D polymer simulation engine:

✅ Contact matrix to restraints conversion
✅ Insulation to backbone stiffness mapping
✅ Langevin dynamics simulation
✅ Reproducibility with fixed seeds
✅ Numerical stability
✅ Chain connectivity constraints
✅ Temperature effects
✅ Edge cases (empty matrices, extreme values)

**Example:**
```bash
pytest tests/test_polymer_sim.py::TestSimulatePolymer -v
```

## Fixtures (conftest.py)

Shared test data and mocks available to all tests:

| Fixture | Description |
|---------|-------------|
| `mock_contact_matrix` | 100x100 synthetic Hi-C matrix with TAD structure |
| `mock_insulation_scores` | Synthetic insulation scores with boundaries |
| `mock_cooler` | Mock cooler.Cooler object for testing |
| `sample_bins_df` | Sample genomic bins DataFrame |
| `sample_insulation_table` | Sample insulation table from cooltools |
| `sample_boundary_df` | Sample boundary DataFrame with prominence |
| `sample_restraints` | Sample polymer restraints for testing |
| `random_seed` | Fixed random seed (42) for reproducibility |

**Using fixtures:**
```python
@pytest.mark.unit
def test_with_mock_cooler(mock_cooler):
    result = your_function(mock_cooler)
    assert result is not None
```

## Writing New Tests

### Test Template

```python
import pytest
from src.your_module import your_function

class TestYourFunction:
    """Test suite for your_function."""

    @pytest.mark.unit
    def test_basic_functionality(self):
        """Test basic functionality works."""
        result = your_function(test_input)
        assert result == expected_output

    @pytest.mark.unit
    def test_edge_case(self):
        """Test edge case handling."""
        with pytest.raises(ValueError):
            your_function(invalid_input)

    @pytest.mark.integration
    @pytest.mark.requires_data
    def test_with_real_data(self):
        """Test with real data file."""
        if not data_available():
            pytest.skip("Test data not found")
        result = your_function(real_data)
        assert result is not None
```

### Best Practices

1. **One test, one behavior** - Each test should verify one specific thing
2. **Descriptive names** - Name tests to describe what they test
3. **Use fixtures** - Leverage conftest.py fixtures for common test data
4. **Mark appropriately** - Use @pytest.mark.unit, .integration, .slow
5. **Test edge cases** - Include tests for empty inputs, NaNs, extreme values
6. **Check coverage** - Aim for >80% coverage on critical modules
7. **Make reproducible** - Use fixed seeds for tests with randomness
8. **Fast by default** - Keep unit tests fast (<1s each)

## Continuous Integration

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: conda-incubator/setup-miniconda@v2
        with:
          environment-file: environment.yml
          activate-environment: hic-analysis
      - name: Run tests
        shell: bash -l {0}
        run: |
          make test-coverage
      - name: Upload coverage
        uses: codecov/codecov-action@v3
```

## Debugging Tests

### Common Issues

**Import errors:**
```bash
# Make sure you're in the project root
cd /path/to/HiC-TAD-Library
pytest
```

**Missing dependencies:**
```bash
conda activate hic-analysis
conda env update -f environment.yml
```

**Tests fail randomly:**
- Check for missing `seed` parameters in random operations
- Use fixtures from conftest.py for consistent test data

### Debug Strategies

```bash
# Run single test with full traceback
pytest tests/test_file.py::test_name -vv --tb=long

# Drop into debugger on failure
pytest tests/test_file.py::test_name --pdb

# Show print statements
pytest tests/test_file.py::test_name -s

# Run last failed tests
pytest --lf

# Run in verbose mode with warnings
pytest -vv -W all
```

## Performance

Test suite performance metrics:

- **Unit tests**: ~0.1-0.5s each
- **Full suite**: ~10-30s (depends on system)
- **Coverage analysis**: +2-5s overhead

Optimize test runs:
```bash
# Skip slow tests during development
pytest -m "not slow"

# Run only tests that changed
pytest --lf

# Run tests in parallel (requires pytest-xdist)
pytest -n auto
```

## CI/CD Integration

The test suite is designed for easy CI/CD integration:

- ✅ No GUI dependencies for core tests
- ✅ Mocked data for most tests (no large downloads)
- ✅ Fast execution (<30s typical)
- ✅ Clear pass/fail signals
- ✅ Coverage reports in XML/HTML format

## Resources

- **[tests/README.md](tests/README.md)** - Detailed test documentation
- **[pytest documentation](https://docs.pytest.org/)** - Official pytest docs
- **[pytest-cov](https://pytest-cov.readthedocs.io/)** - Coverage plugin docs
- **[CLAUDE.md](CLAUDE.md)** - Project architecture and structure
- **[README.md](README.md)** - Project overview

## Support

If tests fail or you encounter issues:

1. Check test output for specific error messages
2. Review [tests/README.md](tests/README.md) for detailed guidance
3. Ensure dependencies are up to date: `conda env update -f environment.yml`
4. Check that you're in the project root directory
5. Try running a single test in verbose mode: `pytest path/to/test -vv`

---

**Quick Reference:**

```bash
make test              # Run all tests
make test-coverage     # Run with coverage report
make test-unit         # Run only unit tests (fast)
pytest -m unit         # Unit tests only
pytest -k "boundary"   # Tests matching "boundary"
pytest --lf            # Run last failed
open htmlcov/index.html # View coverage report
```

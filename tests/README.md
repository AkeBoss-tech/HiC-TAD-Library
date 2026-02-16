# HiC-TAD-Library Test Suite

Comprehensive test suite for the HiC-TAD-Library project.

## Quick Start

```bash
# Run all tests
make test

# Or directly with pytest
pytest

# Run with coverage report
make test-coverage
```

## Test Organization

### Test Files

- **`test_loaders.py`** - Tests for data loading utilities ([src/loaders.py](../src/loaders.py))
  - Loading cooler/mcool files
  - Path construction and resolution handling
  - Error handling for missing files

- **`test_tad_boundaries.py`** - Tests for TAD boundary detection ([src/tad_boundaries.py](../src/tad_boundaries.py))
  - Coordinate parsing and helper functions
  - Directionality Index computation
  - Boundary prominence scoring
  - TAD interval calling
  - Multi-scale insulation analysis
  - Boundary persistence classification
  - Boundary pileup analysis

- **`test_polymer_sim.py`** - Tests for 3D polymer simulation ([src/polymer_sim.py](../src/polymer_sim.py))
  - Contact matrix to restraints conversion
  - Insulation to backbone stiffness mapping
  - Langevin dynamics simulation
  - Numerical stability tests

- **`conftest.py`** - Shared pytest fixtures
  - Mock cooler objects
  - Synthetic contact matrices
  - Sample insulation scores
  - Sample boundary data

## Running Tests

### All Tests
```bash
pytest
# or
make test
```

### Unit Tests Only
```bash
pytest -m unit
# or
make test-unit
```

### Integration Tests Only
```bash
pytest -m integration
# or
make test-integration
```

### Specific Test File
```bash
pytest tests/test_loaders.py
pytest tests/test_tad_boundaries.py
pytest tests/test_polymer_sim.py
```

### Specific Test Class
```bash
pytest tests/test_loaders.py::TestGetCooler
pytest tests/test_tad_boundaries.py::TestBoundaryProminence
```

### Specific Test Function
```bash
pytest tests/test_loaders.py::TestGetCooler::test_get_cooler_default_resolution
```

### Verbose Output
```bash
pytest -vv
# or
make test-verbose
```

## Coverage Reports

### Generate Coverage Report
```bash
make test-coverage
```

This generates:
- Terminal summary with line-by-line coverage
- HTML report in `htmlcov/` directory

### View HTML Coverage Report
```bash
open htmlcov/index.html  # macOS
xdg-open htmlcov/index.html  # Linux
```

## Test Markers

Tests are organized with pytest markers:

- `@pytest.mark.unit` - Fast unit tests that don't require external data
- `@pytest.mark.integration` - Integration tests that may require data files
- `@pytest.mark.slow` - Tests that take longer to run
- `@pytest.mark.requires_data` - Tests that require actual Hi-C data files

### Running Tests by Marker
```bash
# Only fast unit tests
pytest -m unit

# Exclude slow tests
pytest -m "not slow"

# Only tests that require data
pytest -m requires_data

# Combine markers
pytest -m "unit and not slow"
```

## Writing New Tests

### Test Structure Template

```python
import pytest
from src.your_module import your_function

class TestYourFunction:
    """Test suite for your_function."""

    @pytest.mark.unit
    def test_basic_functionality(self):
        """Test that basic functionality works."""
        result = your_function(input_data)
        assert result == expected_output

    @pytest.mark.unit
    def test_edge_case(self):
        """Test edge case handling."""
        result = your_function(edge_case_input)
        assert result is not None
```

### Using Fixtures

```python
@pytest.mark.unit
def test_with_mock_cooler(mock_cooler):
    """Test using the mock cooler fixture."""
    result = your_function(mock_cooler)
    assert result is not None

@pytest.mark.unit
def test_with_contact_matrix(mock_contact_matrix):
    """Test using synthetic contact matrix."""
    result = process_matrix(mock_contact_matrix)
    assert result.shape == (100, 100)
```

## Continuous Integration

To integrate with CI/CD pipelines:

```yaml
# Example GitHub Actions workflow
- name: Run tests
  run: |
    conda activate hic-analysis
    pytest --cov=src --cov-report=xml

- name: Upload coverage
  uses: codecov/codecov-action@v3
```

## Debugging Failed Tests

### Run with Debug Output
```bash
pytest -vv --tb=long
pytest -vv --pdb  # Drop into debugger on failure
```

### Run Specific Failed Test
```bash
pytest tests/test_loaders.py::TestGetCooler::test_specific_test -vv
```

### Print Statements in Tests
```python
def test_debug_example():
    value = calculate_something()
    print(f"Debug: value={value}")  # Will show with pytest -s
    assert value > 0
```

Run with `-s` flag to see print statements:
```bash
pytest -s tests/test_loaders.py
```

## Test Data

Tests use synthetic/mock data generated in `conftest.py`:
- No large Hi-C files required for most tests
- Integration tests marked with `@pytest.mark.requires_data` need actual data
- Skip data-dependent tests if files are missing

## Dependencies

Test dependencies (from [environment.yml](../environment.yml)):
- `pytest` - Test framework
- `pytest-cov` - Coverage plugin
- `pytest-mock` - Mocking utilities

Install with:
```bash
conda env update -f environment.yml
```

## Best Practices

1. **Keep tests fast** - Use mocks instead of real data when possible
2. **Test one thing** - Each test should verify one specific behavior
3. **Use descriptive names** - Test names should describe what they test
4. **Use markers** - Mark tests appropriately (unit, integration, slow)
5. **Check coverage** - Aim for >80% code coverage on critical modules
6. **Test edge cases** - Include tests for empty inputs, NaNs, extreme values
7. **Make tests reproducible** - Use fixed random seeds when needed

## Troubleshooting

### Import Errors
If you get import errors, ensure you're in the project root directory:
```bash
cd /path/to/HiC-TAD-Library
pytest
```

### Missing Dependencies
```bash
conda activate hic-analysis
conda env update -f environment.yml
```

### Permission Errors
Ensure test files are readable:
```bash
chmod -R u+r tests/
```

## Further Reading

- [pytest documentation](https://docs.pytest.org/)
- [pytest-cov documentation](https://pytest-cov.readthedocs.io/)
- [Project README](../README.md)
- [CLAUDE.md](../CLAUDE.md) - Project structure and architecture

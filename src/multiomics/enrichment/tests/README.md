# Enrichment Module Tests

This directory contains comprehensive unit tests for the multiomics enrichment analysis module. The tests cover all major functionality including Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), Gene Ontology (GO) analysis, and Pathway Analysis of Gene Expression (PAGE).

## Test Structure

```
tests/
├── __init__.py              # Test package initialization
├── test_utils.py            # Tests for utility functions
├── test_ora.py              # Tests for Over-Representation Analysis
├── test_gsea.py             # Tests for GSEA functionality
├── test_go_analysis.py      # Tests for GO enrichment analysis
├── test_pypage_wrapper.py   # Tests for PAGE analysis
├── run_tests.py             # Test runner script
└── README.md                # This file
```

## Running Tests

### Quick Start

Run all tests with the test runner:

```bash
# From the tests directory
python run_tests.py

# Or from the enrichment module directory
python -m tests.run_tests
```

### Test Runner Options

The test runner provides several options for different testing scenarios:

```bash
# Run basic tests only (skip integration tests)
python run_tests.py --no-integration

# Include network-dependent tests (GO enrichment)
python run_tests.py --include-network

# Quiet mode (minimal output)
python run_tests.py --quiet

# Verbose mode (detailed output)
python run_tests.py -v -v
```

### Running Individual Test Modules

You can also run individual test modules using unittest:

```bash
# Test utility functions
python -m unittest tests.test_utils

# Test ORA functionality
python -m unittest tests.test_ora

# Test GSEA functionality
python -m unittest tests.test_gsea

# Test GO analysis
python -m unittest tests.test_go_analysis

# Test PAGE analysis
python -m unittest tests.test_pypage_wrapper
```

### Running Specific Test Cases

```bash
# Run specific test class
python -m unittest tests.test_ora.TestORA

# Run specific test method
python -m unittest tests.test_ora.TestORA.test_fisher_exact_test_basic
```

## Test Coverage

### Core Functionality (Always Tested)
- **Utility Functions** (`test_utils.py`)
  - Gene list preparation and cleaning
  - Multiple testing correction (Bonferroni, FDR)
  - Dependency checking
  - Result formatting and validation

### Dependency-Based Testing

The tests automatically detect available dependencies and adjust coverage accordingly:

#### ORA Tests (`test_ora.py`) - Requires `scipy`
- Fisher's exact test
- Hypergeometric test
- Custom ORA implementation
- Multiple testing correction
- Edge cases and error handling

#### GSEA Tests (`test_gsea.py`) - Core + Optional `gseapy`
- **Core functionality** (always tested):
  - Ranked gene list creation
  - Expression matrix preparation
  - Different ranking methods (log2FC, t-test, signal-to-noise)
- **Enhanced functionality** (requires gseapy):
  - Preranked GSEA analysis
  - Full GSEA analysis

#### GO Analysis Tests (`test_go_analysis.py`) - Requires `gseapy` or `goatools`
- **Basic tests** (always run):
  - Method availability checking
  - Gene list validation
- **Network-dependent tests** (optional):
  - GO enrichment with gseapy
  - Comparative analysis between gene sets

#### PAGE Tests (`test_pypage_wrapper.py`) - Requires `pypage`
- Expression profile creation
- Gene set format conversion
- Mutual information analysis
- Permutation testing
- Statistical significance assessment

## Test Categories

### Unit Tests
Test individual functions and methods in isolation:
- Input validation
- Output format verification
- Error handling
- Edge cases

### Integration Tests
Test complete workflows and interactions between components:
- Full analysis pipelines
- Multi-method comparisons
- Realistic biological datasets

### Network Tests
Tests that require internet connectivity (disabled by default):
- GO enrichment via online databases
- Gene set database downloads

## Dependency Management

The test suite gracefully handles missing optional dependencies:

```python
# Tests are automatically skipped when dependencies are unavailable
@unittest.skipUnless(GSEAPY_AVAILABLE, "gseapy not available")
def test_gseapy_functionality(self):
    # This test only runs if gseapy is installed
    pass
```

### Installing Test Dependencies

For complete test coverage, install all optional dependencies:

```bash
# Basic scientific Python stack (required for most tests)
pip install scipy pandas numpy statsmodels

# GSEA and GO analysis
pip install gseapy

# Advanced GO analysis
pip install goatools

# PAGE analysis
pip install bio-pypage
```

## Test Data

Tests use synthetic biological data designed to:
- **Simulate realistic scenarios**: Expression patterns, pathway memberships
- **Provide known outcomes**: Controlled differential expression, pathway enrichment
- **Test edge cases**: Empty inputs, missing genes, invalid parameters
- **Ensure reproducibility**: Fixed random seeds where applicable

### Example Test Data Patterns

```python
# Differential expression simulation
expression_data.loc[upregulated_genes, condition_b_samples] += 2
expression_data.loc[downregulated_genes, condition_a_samples] += 1.5

# Pathway enrichment simulation
gene_sets = {
    'Enriched_Pathway': upregulated_genes[:10],
    'Not_Enriched_Pathway': random_genes,
    'Mixed_Pathway': upregulated_genes[:5] + random_genes[:5]
}
```

## Continuous Integration

The test suite is designed for CI/CD environments:

- **Dependency detection**: Automatically adjusts test coverage
- **Network isolation**: Core tests don't require internet
- **Clear reporting**: Structured output for automated systems
- **Exit codes**: Proper success/failure reporting

### Example CI Usage

```yaml
# GitHub Actions example
- name: Run enrichment tests
  run: |
    cd src/multiomics/enrichment
    python tests/run_tests.py --no-integration --quiet
```

## Troubleshooting

### Common Issues

1. **ImportError for optional dependencies**
   - Expected behavior: Tests are automatically skipped
   - Solution: Install missing dependencies or run with `--no-integration`

2. **Network timeouts during GO tests**
   - Expected behavior: Tests are skipped with network errors
   - Solution: Run without `--include-network` flag

3. **Random test failures**
   - Possible cause: Network-dependent tests or insufficient permutations
   - Solution: Re-run tests or increase permutation numbers in test data

### Debugging Failed Tests

Run tests with maximum verbosity to see detailed error information:

```bash
python run_tests.py -v -v
```

Run individual failing tests:

```bash
python -m unittest tests.test_module.TestClass.test_method -v
```

## Contributing Tests

When adding new functionality to the enrichment module:

1. **Add corresponding tests** in the appropriate test file
2. **Follow naming conventions**: `test_function_name_scenario`
3. **Include docstrings** explaining what is being tested
4. **Test both success and failure cases**
5. **Use appropriate skip decorators** for optional dependencies
6. **Add integration tests** for complex workflows

### Test Template

```python
def test_new_function_basic(self):
    """Test basic functionality of new_function."""
    # Arrange
    input_data = create_test_data()
    
    # Act
    result = new_function(input_data)
    
    # Assert
    self.assertIsInstance(result, expected_type)
    self.assertEqual(len(result), expected_length)
    # ... more assertions
```

## Performance Considerations

Tests are designed to run quickly while maintaining coverage:
- **Small datasets**: Synthetic data with minimal required size
- **Reduced permutations**: Lower iteration counts for statistical tests
- **Efficient assertions**: Focus on correctness rather than exhaustive validation
- **Conditional complexity**: More thorough tests only when dependencies are available

The complete test suite should run in under 2 minutes on typical hardware with all dependencies installed.
"""
Test runner for the enrichment module tests.

This script provides a convenient way to run all tests in the enrichment module
with proper reporting and optional dependency handling.
"""

import unittest
import sys
import os
import warnings

# Add parent directory to path for imports
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

# Import test modules
from tests.test_utils import TestUtils, TestUtilsIntegration
from tests.test_ora import TestORA, TestORAEdgeCases
from tests.test_gsea import TestGSEA, TestGSEAIntegration
from tests.test_go_analysis import TestGOAnalysis, TestGOAnalysisIntegration
from tests.test_pypage_wrapper import TestPypageWrapper, TestPypageIntegration

# Check dependencies
def check_dependencies():
    """Check availability of optional dependencies."""
    dependencies = {}
    
    # Check scipy (required for ORA)
    try:
        import scipy.stats
        dependencies['scipy'] = True
    except ImportError:
        dependencies['scipy'] = False
    
    # Check gseapy (optional for GSEA and GO)
    try:
        import gseapy
        dependencies['gseapy'] = True
    except ImportError:
        dependencies['gseapy'] = False
    
    # Check goatools (optional for GO)
    try:
        import goatools
        dependencies['goatools'] = True
    except ImportError:
        dependencies['goatools'] = False
    
    # Check pypage (optional for PAGE)
    try:
        import pypage
        dependencies['pypage'] = True
    except ImportError:
        dependencies['pypage'] = False
    
    return dependencies


def create_test_suite(include_integration=True, include_network_tests=False):
    """Create test suite with optional components."""
    suite = unittest.TestSuite()
    
    # Always include basic tests
    suite.addTest(unittest.makeSuite(TestUtils))
    suite.addTest(unittest.makeSuite(TestUtilsIntegration))
    
    # Add ORA tests (requires scipy)
    dependencies = check_dependencies()
    if dependencies['scipy']:
        suite.addTest(unittest.makeSuite(TestORA))
        suite.addTest(unittest.makeSuite(TestORAEdgeCases))
    else:
        print("Warning: Skipping ORA tests (scipy not available)")
    
    # Add GSEA tests
    suite.addTest(unittest.makeSuite(TestGSEA))
    if include_integration:
        suite.addTest(unittest.makeSuite(TestGSEAIntegration))
    
    # Add GO tests (may require network access)
    if include_network_tests:
        suite.addTest(unittest.makeSuite(TestGOAnalysis))
        if include_integration:
            suite.addTest(unittest.makeSuite(TestGOAnalysisIntegration))
    else:
        # Add basic GO tests that don't require network
        suite.addTest(unittest.makeSuite(TestGOAnalysis, 'test_get_available_go_libraries'))
        suite.addTest(unittest.makeSuite(TestGOAnalysis, 'test_validate_go_gene_list'))
    
    # Add PAGE tests
    suite.addTest(unittest.makeSuite(TestPypageWrapper))
    if include_integration:
        suite.addTest(unittest.makeSuite(TestPypageIntegration))
    
    return suite


def print_dependency_report():
    """Print a report of available dependencies."""
    print("Enrichment Module Test Runner")
    print("=" * 40)
    
    dependencies = check_dependencies()
    
    print("Dependency Status:")
    for dep, available in dependencies.items():
        status = "✓ Available" if available else "✗ Not Available"
        print(f"  {dep:10}: {status}")
    
    print("\nTest Coverage:")
    print(f"  Core utilities: Always tested")
    print(f"  ORA analysis  : {'✓' if dependencies['scipy'] else '✗'} (requires scipy)")
    print(f"  GSEA analysis : ✓ (basic) / {'✓' if dependencies['gseapy'] else '✗'} (full)")
    print(f"  GO analysis   : {'✓' if dependencies['gseapy'] or dependencies['goatools'] else '✗'}")
    print(f"  PAGE analysis : {'✓' if dependencies['pypage'] else '✗'} (requires pypage)")
    print()


def run_tests(verbosity=2, include_integration=True, include_network_tests=False):
    """Run the test suite."""
    print_dependency_report()
    
    # Suppress warnings for cleaner output
    warnings.filterwarnings('ignore', category=UserWarning)
    warnings.filterwarnings('ignore', category=FutureWarning)
    
    # Create and run test suite
    suite = create_test_suite(include_integration, include_network_tests)
    runner = unittest.TextTestRunner(verbosity=verbosity)
    
    print("Running tests...")
    print("-" * 40)
    
    result = runner.run(suite)
    
    # Print summary
    print("\n" + "=" * 40)
    print("Test Summary:")
    print(f"Tests run: {result.testsRun}")
    print(f"Failures: {len(result.failures)}")
    print(f"Errors: {len(result.errors)}")
    print(f"Skipped: {len(result.skipped) if hasattr(result, 'skipped') else 'N/A'}")
    
    if result.failures:
        print("\nFailures:")
        for test, traceback in result.failures:
            print(f"  - {test}")
    
    if result.errors:
        print("\nErrors:")
        for test, traceback in result.errors:
            print(f"  - {test}")
    
    success = len(result.failures) == 0 and len(result.errors) == 0
    print(f"\nOverall result: {'✓ PASSED' if success else '✗ FAILED'}")
    
    return success


if __name__ == '__main__':
    import argparse
    
    parser = argparse.ArgumentParser(description='Run enrichment module tests')
    parser.add_argument('--no-integration', action='store_true',
                        help='Skip integration tests')
    parser.add_argument('--include-network', action='store_true',
                        help='Include tests that require network access')
    parser.add_argument('--verbose', '-v', action='count', default=2,
                        help='Increase verbosity')
    parser.add_argument('--quiet', '-q', action='store_true',
                        help='Minimal output')
    
    args = parser.parse_args()
    
    verbosity = 0 if args.quiet else args.verbose
    include_integration = not args.no_integration
    
    success = run_tests(
        verbosity=verbosity,
        include_integration=include_integration,
        include_network_tests=args.include_network
    )
    
    sys.exit(0 if success else 1)
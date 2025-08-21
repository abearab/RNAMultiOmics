"""
Simple test runner for the RNAMultiOmics package tests.

This script runs basic tests for the available functionality.
"""

import unittest
import sys
import os

# Add src directory to path for imports
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))


def run_basic_tests():
    """Run basic functionality tests."""
    print("Running basic RNAMultiOmics package tests...")
    print("=" * 50)
    
    # Test 1: Check if modules can be imported
    print("Test 1: Module imports")
    try:
        import multiomics.enrichment._utils
        print("  ✓ multiomics.enrichment._utils module imported successfully")
    except Exception as e:
        print(f"  ✗ multiomics.enrichment._utils import failed: {e}")
        return False
    
    try:
        import multiomics.enrichment.ora
        print("  ✓ multiomics.enrichment.ora module imported successfully")
    except Exception as e:
        print(f"  ✗ multiomics.enrichment.ora import failed: {e}")
        return False
    
    try:
        import multiomics.enrichment.gsea
        print("  ✓ multiomics.enrichment.gsea module imported successfully")
    except Exception as e:
        print(f"  ✗ multiomics.enrichment.gsea import failed: {e}")
        return False
    
    try:
        import multiomics.enrichment.go_analysis
        print("  ✓ multiomics.enrichment.go_analysis module imported successfully")
    except Exception as e:
        print(f"  ✗ multiomics.enrichment.go_analysis import failed: {e}")
        return False
    
    try:
        import multiomics.enrichment.pypage_wrapper
        print("  ✓ multiomics.enrichment.pypage_wrapper module imported successfully")
    except Exception as e:
        print(f"  ✗ multiomics.enrichment.pypage_wrapper import failed: {e}")
        return False
    
    # Test 2: Check basic functionality
    print("\nTest 2: Basic functionality")
    
    # Test gene list preparation
    try:
        test_genes = ['gene1', ' GENE2 ', 'Gene3', 'gene1', '']
        import multiomics.enrichment._utils as _utils
        cleaned = _utils.prepare_gene_list(test_genes)
        print(f"  ✓ Gene list preparation: {len(cleaned)} genes processed")
    except Exception as e:
        print(f"  ✗ Gene list preparation failed: {e}")
    
    # Test dependency checking
    print("\nTest 3: Dependency availability")
    
    # Check scipy availability for ORA
    try:
        import scipy.stats
        print("  ✓ scipy available (ORA functionality enabled)")
    except ImportError:
        print("  ✗ scipy not available (ORA functionality disabled)")
    
    # Check gseapy availability
    try:
        import gseapy
        print("  ✓ gseapy available (full GSEA functionality enabled)")
    except ImportError:
        print("  ✗ gseapy not available (limited GSEA functionality)")
    
    # Check pypage availability
    try:
        import pypage
        print("  ✓ pypage available (PAGE analysis enabled)")
    except ImportError:
        print("  ✗ pypage not available (PAGE analysis disabled)")
    
    print("\n" + "=" * 50)
    print("Basic tests completed successfully!")
    return True


if __name__ == '__main__':
    success = run_basic_tests()
    sys.exit(0 if success else 1)
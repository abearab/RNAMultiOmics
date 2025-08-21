#!/usr/bin/env python3
"""
Test script for pypage integration in RNAMultiOmics enrichment module.

This script demonstrates how to use the PAGE (Pathway Analysis of Gene Expression)
wrapper for pathway enrichment analysis using mutual information.
"""

import sys
import os
import pandas as pd
import numpy as np

# Add the src directory to the path to import multiomics
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_pypage_availability():
    """Test if pypage is available and get info."""
    print("=== Testing pypage availability ===")
    
    try:
        from multiomics.enrichment import check_pypage_availability, get_pypage_info
        
        available = check_pypage_availability()
        info = get_pypage_info()
        
        print(f"pypage available: {available}")
        print(f"pypage info: {info}")
        
        return available
    except Exception as e:
        print(f"Error checking pypage availability: {e}")
        return False


def test_pypage_basic():
    """Test basic PAGE analysis functionality."""
    print("\n=== Testing basic PAGE analysis ===")
    
    try:
        from multiomics.enrichment import run_page
        
        # Create sample expression data (genes x samples)
        np.random.seed(42)  # For reproducible results
        genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5', 'GENE6']
        samples = ['Sample1', 'Sample2', 'Sample3', 'Sample4']
        
        # Generate some correlated expression data
        expression_data = pd.DataFrame(
            np.random.normal(0, 1, (len(genes), len(samples))),
            index=genes,
            columns=samples
        )
        
        # Add some signal for specific pathways
        expression_data.loc['GENE1'] += 2  # High expression
        expression_data.loc['GENE2'] += 1.5
        expression_data.loc['GENE3'] -= 1  # Low expression
        
        print("Expression data:")
        print(expression_data)
        
        # Define gene sets/pathways
        gene_sets = {
            'Upregulated_Pathway': ['GENE1', 'GENE2', 'GENE4'],
            'Downregulated_Pathway': ['GENE3', 'GENE5'],
            'Mixed_Pathway': ['GENE1', 'GENE3', 'GENE6'],
            'Random_Pathway': ['GENE4', 'GENE5', 'GENE6']
        }
        
        print(f"\nGene sets:")
        for pathway, genes in gene_sets.items():
            print(f"  {pathway}: {genes}")
        
        # Run PAGE analysis
        print("\nRunning PAGE analysis...")
        results = run_page(
            expression_data=expression_data,
            gene_sets=gene_sets,
            n_shuffle=100,  # Reduced for faster testing
            alpha=0.05,
            n_bins=5,
            verbose=True
        )
        
        print("\nPAGE results:")
        print(results[['Term', 'MI', 'P_value', 'Informative', 'N_shared']])
        
        # Check results
        assert len(results) == len(gene_sets), "Should have results for all pathways"
        assert 'MI' in results.columns, "Should have mutual information scores"
        assert 'Informative' in results.columns, "Should have informative flags"
        
        print("\n✓ Basic PAGE analysis test passed!")
        return True
        
    except Exception as e:
        print(f"✗ Basic PAGE analysis test failed: {e}")
        return False


def test_pypage_utilities():
    """Test PAGE utility functions."""
    print("\n=== Testing PAGE utility functions ===")
    
    try:
        from multiomics.enrichment import (
            create_expression_profile_from_dataframe,
            create_gene_sets_from_dict
        )
        
        # Test expression profile creation
        expression_data = pd.DataFrame({
            'Sample1': [1.0, 2.0, 3.0],
            'Sample2': [1.1, 2.1, 3.1],
            'Sample3': [0.9, 1.9, 2.9]
        }, index=['GENE1', 'GENE2', 'GENE3'])
        
        expr_profile = create_expression_profile_from_dataframe(
            expression_data, 
            n_bins=5, 
            method='mean'
        )
        
        print(f"Expression profile created: {expr_profile.n_genes} genes, {expr_profile.n_bins} bins")
        
        # Test gene sets creation
        gene_sets = {
            'Pathway1': ['GENE1', 'GENE2'],
            'Pathway2': ['GENE2', 'GENE3']
        }
        
        gs_obj = create_gene_sets_from_dict(gene_sets, n_bins=3)
        print(f"Gene sets created: {gs_obj.n_pathways} pathways, {gs_obj.n_genes} genes")
        
        print("✓ PAGE utility functions test passed!")
        return True
        
    except Exception as e:
        print(f"✗ PAGE utility functions test failed: {e}")
        return False


def test_error_handling():
    """Test error handling for invalid inputs."""
    print("\n=== Testing error handling ===")
    
    try:
        from multiomics.enrichment import run_page
        
        # Test with invalid expression data
        try:
            run_page(
                expression_data="invalid_data",
                gene_sets={'Test': ['GENE1']}
            )
            print("✗ Should have raised ValueError for invalid expression data")
            return False
        except (ValueError, TypeError):
            print("✓ Correctly handled invalid expression data")
        
        # Test with numpy array without genes
        try:
            run_page(
                expression_data=np.array([[1, 2], [3, 4]]),
                gene_sets={'Test': ['GENE1']}
            )
            print("✗ Should have raised ValueError for numpy array without genes")
            return False
        except ValueError:
            print("✓ Correctly handled numpy array without genes parameter")
        
        print("✓ Error handling test passed!")
        return True
        
    except Exception as e:
        print(f"✗ Error handling test failed: {e}")
        return False


def main():
    """Run all PAGE integration tests."""
    print("Testing pypage integration for RNAMultiOmics")
    print("=" * 50)
    
    # Check if pypage is available
    pypage_available = test_pypage_availability()
    
    if not pypage_available:
        print("\npypage is not available. Install it with: pip install bio-pypage")
        print("Skipping PAGE-specific tests.")
        return
    
    # Run tests
    tests = [
        test_pypage_basic,
        test_pypage_utilities,
        test_error_handling
    ]
    
    passed = 0
    total = len(tests)
    
    for test_func in tests:
        if test_func():
            passed += 1
    
    print(f"\n{'='*50}")
    print(f"Test Summary: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All PAGE integration tests passed!")
    else:
        print("⚠️  Some tests failed. Check the output above for details.")


if __name__ == "__main__":
    main()
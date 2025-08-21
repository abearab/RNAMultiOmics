#!/usr/bin/env python3
"""
Simple test for pypage wrapper functionality without full multiomics dependencies.
"""

import sys
import os
import pandas as pd
import numpy as np

# Add the enrichment module to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'multiomics', 'enrichment'))

def test_pypage_wrapper():
    """Test pypage wrapper functionality."""
    print("Testing pypage wrapper...")
    
    try:
        from pypage_wrapper import check_pypage_availability, run_page, get_pypage_info
        
        # Check availability
        available = check_pypage_availability()
        print(f"pypage available: {available}")
        
        if not available:
            print("pypage not available - skipping tests")
            return False
            
        # Get info
        info = get_pypage_info()
        print(f"pypage info: {info}")
        
        # Create test data
        np.random.seed(42)
        genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']
        samples = ['Sample1', 'Sample2', 'Sample3']
        
        # Create expression data with some signal
        expression_data = pd.DataFrame(
            np.random.normal(0, 1, (len(genes), len(samples))),
            index=genes,
            columns=samples
        )
        
        # Add signal to some genes
        expression_data.loc['GENE1'] += 2  # High expression
        expression_data.loc['GENE2'] += 1.5
        expression_data.loc['GENE3'] -= 1  # Low expression
        
        print("Expression data:")
        print(expression_data.round(2))
        
        # Create gene sets
        gene_sets = {
            'Pathway_A': ['GENE1', 'GENE2'],
            'Pathway_B': ['GENE3', 'GENE4'],
            'Pathway_C': ['GENE1', 'GENE3', 'GENE5']
        }
        
        print(f"\nGene sets: {gene_sets}")
        
        # Run PAGE analysis
        print("\nRunning PAGE analysis...")
        results = run_page(
            expression_data=expression_data,
            gene_sets=gene_sets,
            n_shuffle=50,  # Small for testing
            alpha=0.05,
            n_bins=5,
            verbose=True
        )
        
        print("\nResults:")
        print(results)
        
        # Basic validation
        assert len(results) == len(gene_sets)
        assert 'Term' in results.columns
        assert 'MI' in results.columns
        assert 'P_value' in results.columns
        assert 'Informative' in results.columns
        
        print("\n✓ PAGE wrapper test passed!")
        return True
        
    except Exception as e:
        print(f"✗ PAGE wrapper test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    print("Testing pypage wrapper integration")
    print("=" * 40)
    
    success = test_pypage_wrapper()
    
    if success:
        print("\n🎉 All tests passed!")
    else:
        print("\n⚠️ Tests failed!")
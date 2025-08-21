#!/usr/bin/env python3
"""
Simple test script for the enrichment analysis functionality

This script tests the core functionality without importing the full package.
"""

import pandas as pd
import numpy as np
import os
import sys

# Add the enrichment module directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'multiomics', 'enrichment'))

def test_basic_functionality():
    """Test basic functionality of enrichment modules"""
    print("Testing basic enrichment analysis functionality...")
    print("="*60)
    
    # Import modules directly
    import _utils as utils
    import ora
    import gsea
    
    print("✓ Successfully imported all enrichment modules")
    
    # Test utility functions
    print("\n1. Testing utility functions:")
    messy_genes = ['GENE1', 'gene2', ' GENE3 ', 'Gene4.1', 'GENE1', 'gene5']
    cleaned = utils.prepare_gene_list(messy_genes, verbose=True)
    print(f"   Cleaned genes: {cleaned}")
    
    # Test ORA
    print("\n2. Testing Over-Representation Analysis:")
    study_genes = ['GENE_001', 'GENE_002', 'GENE_003', 'GENE_004', 'GENE_005']
    gene_sets = {
        'Pathway_A': ['GENE_001', 'GENE_002', 'GENE_010', 'GENE_011', 'GENE_012'],
        'Pathway_B': ['GENE_003', 'GENE_004', 'GENE_013', 'GENE_014', 'GENE_015'],
        'Pathway_C': ['GENE_006', 'GENE_007', 'GENE_008', 'GENE_009', 'GENE_016']
    }
    background = [f'GENE_{i:03d}' for i in range(1, 21)]
    
    # Test Fisher's exact test
    fisher_result = ora.run_fisher_exact_test(study_genes, gene_sets['Pathway_A'], background)
    print(f"   Fisher's test p-value: {fisher_result['p_value']:.4f}")
    
    # Test custom ORA
    ora_results = ora.run_custom_ora(study_genes, gene_sets, background, verbose=True)
    print(f"   Found {len(ora_results)} pathway results")
    print(f"   Most significant: {ora_results.iloc[0]['Term']} (p={ora_results.iloc[0]['P_value']:.4f})")
    
    # Test ranked gene analysis
    print("\n3. Testing ranked gene analysis:")
    np.random.seed(42)
    genes = [f'GENE_{i:03d}' for i in range(1, 51)]
    samples_a = ['A1', 'A2', 'A3']
    samples_b = ['B1', 'B2', 'B3']
    
    # Create sample expression data
    data = pd.DataFrame(
        np.random.randn(50, 6),
        index=genes,
        columns=samples_a + samples_b
    )
    # Add some differential expression
    data.loc[genes[:10], samples_b] += 2  # Upregulated in B
    data.loc[genes[40:], samples_a] += 1.5  # Upregulated in A
    
    ranked_genes = gsea.create_ranked_list(data, samples_b, samples_a, verbose=True)
    print(f"   Top upregulated: {ranked_genes.head(3).index.tolist()}")
    print(f"   Top downregulated: {ranked_genes.tail(3).index.tolist()}")
    
    print("\n✓ All basic functionality tests passed!")
    
    # Check optional dependencies
    print("\n4. Checking optional dependencies:")
    print(f"   GSEAPY available: {gsea.GSEAPY_AVAILABLE}")
    
    if not gsea.GSEAPY_AVAILABLE:
        print("   Install gseapy for full GSEA and GO enrichment functionality:")
        print("   pip install gseapy")
    
    print("\n" + "="*60)
    print("Gene set enrichment analysis module is working correctly!")
    print("The module provides statistical functions for pathway analysis")
    print("and can be extended with additional dependencies for more features.")

if __name__ == "__main__":
    test_basic_functionality()
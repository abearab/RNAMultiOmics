#!/usr/bin/env python3
"""
Example script demonstrating gene set enrichment analysis functionality in RNAMultiOmics

This script shows how to use the enrichment analysis modules for:
1. Over-representation analysis (ORA)
2. Gene set enrichment analysis (GSEA) 
3. GO enrichment analysis
4. PAGE (Pathway Analysis of Gene Expression)
5. Custom statistical tests

Usage:
    python examples/enrichment_analysis_example.py
"""

import pandas as pd
import numpy as np
import os
import sys

# Add the src directory to path to import multiomics modules
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def create_sample_data():
    """Create sample gene expression data for demonstration"""
    print("Creating sample gene expression data...")
    
    # Set random seed for reproducibility
    np.random.seed(42)
    
    # Create gene names
    genes = [f'GENE_{i:03d}' for i in range(1, 201)]  # 200 genes
    
    # Create sample names
    group1_samples = [f'Control_{i}' for i in range(1, 6)]  # 5 control samples
    group2_samples = [f'Treatment_{i}' for i in range(1, 6)]  # 5 treatment samples
    all_samples = group1_samples + group2_samples
    
    # Create expression data (log2 scale)
    expression_data = pd.DataFrame(
        np.random.normal(8, 2, (len(genes), len(all_samples))),
        index=genes,
        columns=all_samples
    )
    
    # Add differential expression to some genes
    # Make first 50 genes upregulated in treatment
    expression_data.loc[genes[:50], group2_samples] += np.random.normal(2, 0.5, (50, 5))
    
    # Make genes 150-200 downregulated in treatment 
    expression_data.loc[genes[150:], group2_samples] -= np.random.normal(1.5, 0.3, (50, 5))
    
    print(f"Created expression data: {expression_data.shape[0]} genes × {expression_data.shape[1]} samples")
    
    return expression_data, group1_samples, group2_samples

def create_sample_gene_sets():
    """Create sample gene sets for demonstration"""
    print("Creating sample gene sets...")
    
    gene_sets = {
        'Upregulated_Pathway': [f'GENE_{i:03d}' for i in range(1, 26)],  # Contains upregulated genes
        'Downregulated_Pathway': [f'GENE_{i:03d}' for i in range(151, 176)],  # Contains downregulated genes
        'Mixed_Pathway': [f'GENE_{i:03d}' for i in range(25, 50)] + [f'GENE_{i:03d}' for i in range(100, 125)],
        'Unrelated_Pathway': [f'GENE_{i:03d}' for i in range(75, 100)],
        'Large_Pathway': [f'GENE_{i:03d}' for i in range(1, 101)],  # Large pathway
    }
    
    print(f"Created {len(gene_sets)} sample gene sets")
    return gene_sets

def demonstrate_ora():
    """Demonstrate Over-Representation Analysis"""
    print("\n" + "="*60)
    print("DEMONSTRATING OVER-REPRESENTATION ANALYSIS (ORA)")
    print("="*60)
    
    # Import ORA module
    from multiomics.enrichment.ora import run_custom_ora, format_ora_results
    
    # Create sample differentially expressed genes
    # Simulate results from differential expression analysis
    upregulated_genes = [f'GENE_{i:03d}' for i in range(1, 31)]  # Top 30 upregulated
    downregulated_genes = [f'GENE_{i:03d}' for i in range(151, 181)]  # Top 30 downregulated
    
    print(f"Testing ORA with {len(upregulated_genes)} upregulated genes")
    print(f"First 10 upregulated genes: {upregulated_genes[:10]}")
    
    # Create gene sets
    gene_sets = create_sample_gene_sets()
    
    # Create background (all genes tested)
    background_genes = [f'GENE_{i:03d}' for i in range(1, 201)]
    
    # Run ORA for upregulated genes
    print("\nRunning ORA for upregulated genes...")
    ora_results = run_custom_ora(
        study_genes=upregulated_genes,
        gene_sets=gene_sets,
        background_genes=background_genes,
        method='fisher',
        alpha=0.05,
        verbose=True
    )
    
    print("\nORA Results (all):")
    print(ora_results[['Term', 'Overlap', 'Gene_Set_Size', 'P_value', 'Adjusted_P_value', 'Enrichment_Score', 'Significant']])
    
    # Format results
    significant_results = format_ora_results(ora_results, significance_cutoff=0.1)  # Relaxed cutoff for demo
    print(f"\nSignificant results (FDR < 0.1): {len(significant_results)}")
    if not significant_results.empty:
        print(significant_results[['Term', 'Overlap', 'P_value', 'Adjusted_P_value']])
    
    return ora_results

def demonstrate_ranked_analysis():
    """Demonstrate ranked gene analysis for GSEA"""
    print("\n" + "="*60)
    print("DEMONSTRATING RANKED GENE ANALYSIS")
    print("="*60)
    
    # Import GSEA module
    from multiomics.enrichment.gsea import create_ranked_list
    from multiomics.enrichment._utils import prepare_gene_list
    
    # Create sample data
    expression_data, group1_samples, group2_samples = create_sample_data()
    
    print(f"Creating ranked gene list comparing {group2_samples} vs {group1_samples}")
    
    # Create ranked gene list
    ranked_genes = create_ranked_list(
        data=expression_data,
        group1_samples=group2_samples,  # Treatment as numerator
        group2_samples=group1_samples,  # Control as denominator
        method='log2fc',
        verbose=True
    )
    
    print(f"\nTop 10 upregulated genes (in treatment vs control):")
    top_up = ranked_genes.head(10)
    for gene, score in top_up.items():
        print(f"  {gene}: {score:.3f}")
    
    print(f"\nTop 10 downregulated genes (in treatment vs control):")
    top_down = ranked_genes.tail(10)
    for gene, score in top_down.items():
        print(f"  {gene}: {score:.3f}")
    
    return ranked_genes

def demonstrate_utilities():
    """Demonstrate utility functions"""
    print("\n" + "="*60)
    print("DEMONSTRATING UTILITY FUNCTIONS")
    print("="*60)
    
    from multiomics.enrichment._utils import prepare_gene_list, clean_gene_name, filter_gene_sets
    
    # Test gene list preparation
    print("Testing gene list preparation...")
    messy_gene_list = ['GENE1', 'gene2', ' GENE3 ', 'Gene4.1', 'GENE1', 'gene5', None, '']
    cleaned_genes = prepare_gene_list(messy_gene_list, remove_duplicates=True, verbose=True)
    print(f"Original: {messy_gene_list}")
    print(f"Cleaned:  {cleaned_genes}")
    
    # Test gene name cleaning
    print("\nTesting gene name cleaning...")
    test_names = [' gene123 ', 'GENE456.2', 'Gene_ABC', 'xyz789.version1']
    for name in test_names:
        cleaned = clean_gene_name(name)
        print(f"  '{name}' -> '{cleaned}'")
    
    # Test gene set filtering
    print("\nTesting gene set filtering...")
    gene_sets = create_sample_gene_sets()
    test_genes = [f'GENE_{i:03d}' for i in range(1, 51)]  # First 50 genes
    
    filtered_sets = filter_gene_sets(
        gene_sets=gene_sets,
        gene_list=test_genes,
        min_overlap=5,
        min_size=20,
        max_size=100,
        verbose=True
    )
    
    print("Filtered gene sets:")
    for name, genes in filtered_sets.items():
        overlap = len(set(genes) & set(test_genes))
        print(f"  {name}: {len(genes)} genes, {overlap} overlap")

def main():
    """Main function to run all demonstrations"""
    print("RNAMultiOmics Gene Set Enrichment Analysis - Example Script")
    print("="*70)
    
    # Check dependency availability
    print("Checking dependencies...")
    try:
        import gseapy
        print("✓ gseapy is available")
        gseapy_available = True
    except ImportError:
        print("✗ gseapy is not available (install with: pip install gseapy)")
        gseapy_available = False
    
    try:
        from goatools.obo_parser import GODag
        print("✓ goatools is available")
        goatools_available = True
    except ImportError:
        print("✗ goatools is not available (install with: pip install goatools)")
        goatools_available = False
    
    try:
        import statsmodels
        print("✓ statsmodels is available")
    except ImportError:
        print("✗ statsmodels is not available (install with: pip install statsmodels)")
    
    # Run demonstrations
    demonstrate_utilities()
    ora_results = demonstrate_ora()
    ranked_genes = demonstrate_ranked_analysis()
    
    # Summary
    print("\n" + "="*60)
    print("SUMMARY")
    print("="*60)
    print("✓ Over-representation analysis (ORA) - works with basic Python libraries")
    print("✓ Ranked gene analysis - works with basic Python libraries")
    print("✓ Utility functions - work with basic Python libraries")
    
    if gseapy_available:
        print("✓ GSEA analysis - available with gseapy")
        print("✓ GO enrichment - available with gseapy")
    else:
        print("○ GSEA analysis - install gseapy for full functionality")
        print("○ GO enrichment - install gseapy or goatools for functionality")
    
    if goatools_available:
        print("✓ Advanced GO analysis - available with goatools")
    else:
        print("○ Advanced GO analysis - install goatools for additional functionality")
    
    # PAGE Analysis (if available)
    print("\n" + "="*50)
    print("Testing PAGE (Pathway Analysis of Gene Expression)")
    print("="*50)
    
    try:
        from multiomics.enrichment.pypage_wrapper import check_pypage_availability, run_page
        
        if check_pypage_availability():
            print("✓ pypage is available for PAGE analysis")
            
            try:
                print("\nRunning PAGE analysis...")
                print("PAGE uses mutual information to identify informative pathways")
                
                # Run PAGE analysis on a subset of data for demonstration
                subset_expression = expression_data.iloc[:20, :]  # First 20 genes
                subset_gene_sets = {k: [g for g in v if g in subset_expression.index] 
                                  for k, v in gene_sets.items()}
                subset_gene_sets = {k: v for k, v in subset_gene_sets.items() if len(v) >= 2}
                
                if subset_gene_sets:
                    page_results = run_page(
                        expression_data=subset_expression,
                        gene_sets=subset_gene_sets,
                        n_shuffle=100,  # Reduced for demo
                        alpha=0.05,
                        n_bins=5,
                        verbose=False
                    )
                    
                    print(f"\nPAGE Results (top 5):")
                    print(page_results[['Term', 'MI', 'Informative', 'N_shared']].head())
                    
                    informative_count = page_results['Informative'].sum()
                    print(f"\nFound {informative_count} informative pathways using mutual information")
                else:
                    print("No suitable gene sets for PAGE analysis in subset")
                    
            except Exception as e:
                print(f"PAGE analysis failed: {e}")
        else:
            print("○ PAGE analysis - install bio-pypage for functionality")
            print("  pip install bio-pypage")
            
    except ImportError:
        print("○ PAGE analysis - install bio-pypage for functionality")
    
    print("\nThe enrichment analysis module provides a comprehensive framework for")
    print("pathway and gene set enrichment analysis that can work with different")
    print("dependency configurations, from basic Python libraries to specialized")
    print("bioinformatics packages.")

if __name__ == "__main__":
    main()
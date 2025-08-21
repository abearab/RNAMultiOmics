#!/usr/bin/env python3
"""
Integration test to verify the enrichment module works with sample expression data
"""

import pandas as pd
import numpy as np
import sys
import os

# Add src to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

def test_with_expression_data():
    """Test enrichment analysis with realistic expression data"""
    print("Testing enrichment analysis integration...")
    
    # Create realistic sample data similar to what would come from expression analysis
    np.random.seed(42)
    
    # Gene names (using common gene symbols)
    genes = [
        'GAPDH', 'ACTB', 'TP53', 'MYC', 'EGFR', 'KRAS', 'PIK3CA', 'AKT1', 'MTOR', 'PTEN',
        'BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'MLH1', 'MSH2', 'MSH6', 'PMS2', 'APC',
        'CDKN2A', 'RB1', 'CCND1', 'CDK4', 'CDK6', 'E2F1', 'E2F3', 'MYB', 'PCNA', 'CCNE1',
        'TGFB1', 'SMAD2', 'SMAD3', 'SMAD4', 'TGFBR1', 'TGFBR2', 'BMP2', 'BMP4', 'NOGGIN', 'GREMLIN1',
        'WNT1', 'WNT3A', 'CTNNB1', 'APC2', 'AXIN1', 'AXIN2', 'GSK3B', 'DVL1', 'DVL2', 'DVL3'
    ]
    
    # Sample conditions
    control_samples = ['Control_1', 'Control_2', 'Control_3', 'Control_4']
    treatment_samples = ['Treatment_1', 'Treatment_2', 'Treatment_3', 'Treatment_4']
    
    # Create expression matrix (log2 TPM values)
    expression_data = pd.DataFrame(
        np.random.normal(5, 2, (len(genes), len(control_samples + treatment_samples))),
        index=genes,
        columns=control_samples + treatment_samples
    )
    
    # Add differential expression patterns
    # Upregulate cancer-related genes in treatment
    cancer_genes = ['TP53', 'MYC', 'EGFR', 'KRAS', 'PIK3CA', 'AKT1', 'MTOR']
    expression_data.loc[cancer_genes, treatment_samples] += np.random.normal(2, 0.5, (len(cancer_genes), len(treatment_samples)))
    
    # Downregulate tumor suppressors
    tumor_suppressors = ['BRCA1', 'BRCA2', 'ATM', 'PTEN', 'RB1', 'CDKN2A']
    expression_data.loc[tumor_suppressors, treatment_samples] -= np.random.normal(1.5, 0.3, (len(tumor_suppressors), len(treatment_samples)))
    
    print(f"Created expression data: {expression_data.shape[0]} genes × {expression_data.shape[1]} samples")
    
    # Import enrichment modules directly
    sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src', 'multiomics', 'enrichment'))
    import gsea, ora, _utils
    
    print("\n1. Creating ranked gene list...")
    ranked_genes = gsea.create_ranked_list(
        data=expression_data,
        group1=treatment_samples,
        group2=control_samples,
        method='log2fc',
        verbose=True
    )
    
    print(f"Top 5 upregulated genes: {ranked_genes.head().index.tolist()}")
    print(f"Top 5 downregulated genes: {ranked_genes.tail().index.tolist()}")
    
    print("\n2. Testing ORA with pathway gene sets...")
    
    # Define biologically relevant gene sets
    pathway_gene_sets = {
        'p53_Pathway': ['TP53', 'ATM', 'CHEK2', 'CDKN2A', 'RB1', 'MDM2', 'MDM4'],
        'PI3K_AKT_Pathway': ['PIK3CA', 'AKT1', 'MTOR', 'PTEN', 'TSC1', 'TSC2'],
        'Cell_Cycle': ['CCND1', 'CDK4', 'CDK6', 'E2F1', 'E2F3', 'PCNA', 'CCNE1', 'RB1'],
        'DNA_Repair': ['BRCA1', 'BRCA2', 'ATM', 'CHEK2', 'PALB2', 'MLH1', 'MSH2', 'MSH6'],
        'Housekeeping': ['GAPDH', 'ACTB', 'TUBB', 'HPRT1', 'GUSB']
    }
    
    # Get significantly upregulated genes (top 20% with positive log2FC)
    upregulated_threshold = ranked_genes.quantile(0.8)
    upregulated_genes = ranked_genes[ranked_genes > upregulated_threshold].index.tolist()
    
    print(f"Upregulated genes (top 20%): {upregulated_genes}")
    
    # Run ORA
    ora_results = ora.run_custom_ora(
        study_genes=upregulated_genes,
        gene_sets=pathway_gene_sets,
        background_genes=genes,
        method='fisher',
        alpha=0.05,
        verbose=True
    )
    
    print("\nORA Results:")
    print(ora_results[['Term', 'Overlap', 'Gene_Set_Size', 'P_value', 'Adjusted_P_value', 'Enrichment_Score']])
    
    print("\n3. Testing utility functions...")
    
    # Test gene list preparation
    messy_genes = ['tp53', ' BRCA1 ', 'brca2.1', 'AKT1', 'akt1']
    cleaned = _utils.prepare_gene_list(messy_genes, remove_duplicates=True, verbose=True)
    print(f"Cleaned gene list: {cleaned}")
    
    print("\n✓ Integration test completed successfully!")
    print("\nThe enrichment module successfully:")
    print("- Created ranked gene lists from expression data")
    print("- Performed statistical over-representation analysis") 
    print("- Integrated with realistic biological pathways")
    print("- Handled gene name normalization and cleaning")
    
    return True

if __name__ == "__main__":
    test_with_expression_data()
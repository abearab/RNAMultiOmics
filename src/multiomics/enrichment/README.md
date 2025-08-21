# Gene Set Enrichment Analysis Module

This module provides comprehensive gene set enrichment analysis functionality for the RNAMultiOmics package, including Over-Representation Analysis (ORA), Gene Set Enrichment Analysis (GSEA), and Gene Ontology (GO) enrichment analysis.

## Features

### Core Statistical Functions
- **Over-Representation Analysis (ORA)**: Fisher's exact test and hypergeometric test for gene set enrichment
- **Ranked Gene Analysis**: Create ranked gene lists from expression data for GSEA
- **Statistical Testing**: Multiple testing correction using FDR (Benjamini-Hochberg)
- **Utility Functions**: Gene list preparation, cleaning, and formatting

### Integration with Popular Tools
- **gseapy**: GSEA, pre-ranked GSEA, and GO enrichment analysis
- **goatools**: Advanced GO enrichment analysis with custom gene associations
- **Custom implementations**: Statistical tests that work with basic Python libraries

## Installation

The enrichment module works with basic Python scientific libraries but can be enhanced with additional packages:

### Required Dependencies (already included in RNAMultiOmics)
```bash
pip install pandas numpy scipy statsmodels
```

### Optional Dependencies for Enhanced Functionality
```bash
# For GSEA and online GO enrichment
pip install gseapy

# For advanced GO analysis with custom databases
pip install goatools

# For extended bioinformatics functionality
pip install bioservices
```

## Quick Start

### Basic Over-Representation Analysis

```python
from multiomics.enrichment import run_ora, run_custom_ora

# Your differentially expressed genes
de_genes = ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE5']

# Define gene sets (pathways)
gene_sets = {
    'Pathway_A': ['GENE1', 'GENE2', 'GENE10', 'GENE11'],
    'Pathway_B': ['GENE3', 'GENE4', 'GENE12', 'GENE13'],
    'Pathway_C': ['GENE6', 'GENE7', 'GENE8', 'GENE9']
}

# Background genes (all genes tested)
background = [f'GENE{i}' for i in range(1, 101)]

# Run ORA
results = run_custom_ora(
    study_genes=de_genes,
    gene_sets=gene_sets, 
    background_genes=background,
    method='fisher',
    verbose=True
)

print(results[['Term', 'P_value', 'Adjusted_P_value', 'Significant']])
```

### Gene Set Enrichment Analysis (GSEA)

```python
from multiomics.enrichment import run_gsea, run_prerank, create_ranked_list

# Create ranked gene list from expression data
ranked_genes = create_ranked_list(
    data=expression_df,  # genes x samples DataFrame
    group1_samples=['Treatment1', 'Treatment2', 'Treatment3'],
    group2_samples=['Control1', 'Control2', 'Control3'],
    method='log2fc'
)

# Run pre-ranked GSEA (requires gseapy)
if gseapy_available:
    gsea_results = run_prerank(
        rnk=ranked_genes,
        gene_sets='KEGG_2016',  # or custom gene sets
        permutation_num=1000
    )
```

### GO Enrichment Analysis

```python
from multiomics.enrichment import run_go_enrichment

# Run GO enrichment (requires gseapy or goatools)
go_results = run_go_enrichment(
    gene_list=de_genes,
    method='gseapy',  # or 'goatools'
    organism='Human',
    cutoff=0.05
)
```

### Utility Functions

```python
from multiomics.enrichment import prepare_gene_list, format_results

# Clean and prepare gene lists
messy_genes = ['gene1', ' GENE2 ', 'Gene3.1', 'GENE1']  # duplicates, case issues
clean_genes = prepare_gene_list(messy_genes, remove_duplicates=True)

# Format results for output
formatted = format_results(results, top_n=20, significance_cutoff=0.05)
```

## Module Structure

```
src/multiomics/enrichment/
├── __init__.py          # Main module interface
├── gsea.py             # GSEA and ranked analysis functions
├── go_analysis.py      # GO enrichment analysis
├── ora.py              # Over-representation analysis
└── _utils.py           # Utility functions
```

## Function Reference

### ORA Functions
- `run_ora()`: Over-representation analysis using gseapy
- `run_custom_ora()`: Custom ORA with Fisher's exact or hypergeometric test
- `run_fisher_exact_test()`: Individual Fisher's exact test
- `run_hypergeometric_test()`: Individual hypergeometric test

### GSEA Functions  
- `run_gsea()`: Full GSEA analysis
- `run_prerank()`: Pre-ranked GSEA
- `create_ranked_list()`: Create ranked gene lists from expression data

### GO Analysis Functions
- `run_go_enrichment()`: GO enrichment with multiple backends
- `run_go_enrichment_gseapy()`: GO enrichment using gseapy
- `run_go_enrichment_goatools()`: GO enrichment using goatools

### Utility Functions
- `prepare_gene_list()`: Clean and prepare gene lists
- `clean_gene_name()`: Standardize gene names
- `filter_gene_sets()`: Filter gene sets by size and overlap
- `format_results()`: Standardize result formatting
- `extract_genes_from_anndata()`: Extract data from AnnData objects

## Examples

See the `examples/` directory for complete examples:
- `test_enrichment_basic.py`: Basic functionality test
- `enrichment_analysis_example.py`: Comprehensive example (requires full dependencies)

## Dependency Management

The module is designed to work gracefully with different dependency configurations:

1. **Basic Mode**: Works with pandas, numpy, scipy, statsmodels
   - Custom statistical tests
   - Gene list utilities
   - Basic enrichment analysis

2. **Enhanced Mode**: Add gseapy for full functionality
   - Online gene set databases
   - GSEA analysis
   - GO enrichment

3. **Advanced Mode**: Add goatools for specialized GO analysis
   - Custom GO databases
   - Advanced GO statistics
   - Local GO analysis

## Gene Set Databases

When using gseapy, you can access many online databases:
- `'KEGG_2016'`, `'KEGG_2021'`: KEGG pathway databases
- `'GO_Biological_Process_2023'`: GO biological processes
- `'GO_Cellular_Component_2023'`: GO cellular components  
- `'GO_Molecular_Function_2023'`: GO molecular functions
- `'Reactome_2022'`: Reactome pathway database
- `'WikiPathways_2019'`: WikiPathways database

## Tips and Best Practices

1. **Gene ID Consistency**: Ensure gene IDs match between your data and gene sets
2. **Background Selection**: Use appropriate background gene sets (all tested genes)
3. **Multiple Testing**: Always apply multiple testing correction for pathway analysis
4. **Sample Size**: Ensure adequate sample sizes for statistical testing
5. **Gene Set Size**: Filter gene sets to appropriate sizes (typically 15-500 genes)

## Citation

If you use this module in your research, please cite:
- The RNAMultiOmics package
- gseapy if used: Fang et al. (2016) Bioinformatics
- goatools if used: Klopfenstein et al. (2018) Scientific Reports

## Contributing

Contributions are welcome! Please see the main repository for contribution guidelines.
"""
Utility functions for gene set enrichment analysis
"""

import pandas as pd
import numpy as np
from typing import Union, Optional, List, Dict, Any, Tuple
import re


def prepare_gene_list(
    data: Union[pd.DataFrame, pd.Series, List[str]],
    gene_column: Optional[str] = None,
    value_column: Optional[str] = None,
    convert_to_symbol: bool = False,
    remove_duplicates: bool = True,
    verbose: bool = False
) -> List[str]:
    """
    Prepare gene list for enrichment analysis
    
    Parameters
    ----------
    data : pd.DataFrame, pd.Series, or list
        Input data containing gene identifiers
    gene_column : str, optional
        Column name containing gene identifiers (for DataFrame input)
    value_column : str, optional
        Column name containing values for ranking (for DataFrame input)
    convert_to_symbol : bool, default False
        Whether to attempt conversion to gene symbols
    remove_duplicates : bool, default True
        Whether to remove duplicate gene names
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    list
        Cleaned list of gene identifiers
    """
    if isinstance(data, list):
        gene_list = data.copy()
    elif isinstance(data, pd.Series):
        gene_list = data.tolist()
    elif isinstance(data, pd.DataFrame):
        if gene_column is None:
            # Use index as gene identifiers
            gene_list = data.index.tolist()
        else:
            gene_list = data[gene_column].tolist()
    else:
        raise ValueError("Input data must be a list, pandas Series, or pandas DataFrame")
    
    # Convert to strings and remove any NaN values
    gene_list = [str(gene) for gene in gene_list if pd.notna(gene)]
    
    # Remove duplicates if requested
    if remove_duplicates:
        original_length = len(gene_list)
        gene_list = list(dict.fromkeys(gene_list))  # Preserves order
        if verbose and len(gene_list) != original_length:
            print(f"Removed {original_length - len(gene_list)} duplicate genes")
    
    # Clean gene names (remove extra whitespace, standardize format)
    gene_list = [clean_gene_name(gene) for gene in gene_list]
    
    if verbose:
        print(f"Prepared gene list with {len(gene_list)} genes")
    
    return gene_list


def clean_gene_name(gene_name: str) -> str:
    """
    Clean and standardize gene names
    
    Parameters
    ----------
    gene_name : str
        Raw gene name
        
    Returns
    -------
    str
        Cleaned gene name
    """
    # Remove extra whitespace
    gene_name = str(gene_name).strip()
    
    # Convert to uppercase for gene symbols
    gene_name = gene_name.upper()
    
    # Remove common prefixes/suffixes that might interfere
    # Remove version numbers (e.g., ENSG00000123456.1 -> ENSG00000123456)
    gene_name = re.sub(r'\.\d+$', '', gene_name)
    
    return gene_name


def create_ranked_gene_list(
    data: pd.DataFrame,
    group1_samples: List[str],
    group2_samples: List[str],
    method: str = 'log2fc',
    gene_column: Optional[str] = None,
    verbose: bool = False
) -> pd.Series:
    """
    Create ranked gene list for GSEA from expression data
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data with genes as rows and samples as columns
    group1_samples : list
        Sample names for group 1 (numerator)
    group2_samples : list
        Sample names for group 2 (denominator)
    method : str, default 'log2fc'
        Ranking method: 'log2fc', 'fold_change', 't_test', 'signal2noise'
    gene_column : str, optional
        Column name containing gene identifiers (if genes not in index)
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    pd.Series
        Ranked gene list with gene names as index and ranking values
    """
    from scipy import stats
    
    # Ensure we have the samples in the data
    missing_samples = set(group1_samples + group2_samples) - set(data.columns)
    if missing_samples:
        raise ValueError(f"Missing samples in data: {missing_samples}")
    
    # Get expression data for each group
    group1_data = data[group1_samples]
    group2_data = data[group2_samples]
    
    if method == 'log2fc':
        # Calculate log2 fold change
        mean1 = group1_data.mean(axis=1)
        mean2 = group2_data.mean(axis=1)
        # Add small pseudocount to avoid log(0)
        pseudocount = 1e-8
        ranking_values = np.log2((mean1 + pseudocount) / (mean2 + pseudocount))
        
    elif method == 'fold_change':
        # Calculate fold change
        mean1 = group1_data.mean(axis=1)
        mean2 = group2_data.mean(axis=1)
        ranking_values = mean1 / (mean2 + 1e-8)
        
    elif method == 't_test':
        # Perform t-test for each gene
        t_stats = []
        for gene in data.index:
            t_stat, _ = stats.ttest_ind(group1_data.loc[gene], group2_data.loc[gene])
            t_stats.append(t_stat)
        ranking_values = pd.Series(t_stats, index=data.index)
        
    elif method == 'signal2noise':
        # Calculate signal-to-noise ratio
        mean1 = group1_data.mean(axis=1)
        mean2 = group2_data.mean(axis=1)
        std1 = group1_data.std(axis=1)
        std2 = group2_data.std(axis=1)
        ranking_values = (mean1 - mean2) / (std1 + std2 + 1e-8)
        
    else:
        raise ValueError(f"Unknown ranking method: {method}")
    
    # Sort by ranking values in descending order
    ranked_genes = ranking_values.sort_values(ascending=False)
    
    # Remove any NaN values
    ranked_genes = ranked_genes.dropna()
    
    if verbose:
        print(f"Created ranked gene list with {len(ranked_genes)} genes using {method}")
        print(f"Top 5 upregulated: {ranked_genes.head().index.tolist()}")
        print(f"Top 5 downregulated: {ranked_genes.tail().index.tolist()}")
    
    return ranked_genes


def filter_gene_sets(
    gene_sets: Dict[str, List[str]],
    gene_list: List[str],
    min_overlap: int = 5,
    min_size: int = 15,
    max_size: int = 500,
    verbose: bool = False
) -> Dict[str, List[str]]:
    """
    Filter gene sets based on size and overlap with input gene list
    
    Parameters
    ----------
    gene_sets : dict
        Dictionary of gene sets
    gene_list : list
        List of input genes
    min_overlap : int, default 5
        Minimum overlap with input gene list
    min_size : int, default 15
        Minimum gene set size
    max_size : int, default 500
        Maximum gene set size
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    dict
        Filtered gene sets
    """
    gene_set = set(gene_list)
    filtered_sets = {}
    
    for name, genes in gene_sets.items():
        gene_set_genes = set(genes)
        overlap = len(gene_set & gene_set_genes)
        size = len(gene_set_genes)
        
        if (overlap >= min_overlap and 
            min_size <= size <= max_size):
            filtered_sets[name] = genes
    
    if verbose:
        print(f"Filtered gene sets: {len(filtered_sets)} / {len(gene_sets)} retained")
        print(f"Criteria: min_overlap={min_overlap}, min_size={min_size}, max_size={max_size}")
    
    return filtered_sets


def format_results(
    results: Union[pd.DataFrame, List, Dict],
    source: str = 'gseapy',
    sort_by: str = 'Adjusted P-value',
    ascending: bool = True,
    top_n: Optional[int] = None,
    significance_cutoff: float = 0.05,
    verbose: bool = False
) -> pd.DataFrame:
    """
    Format enrichment analysis results into standardized DataFrame
    
    Parameters
    ----------
    results : pd.DataFrame, list, or dict
        Enrichment analysis results
    source : str, default 'gseapy'
        Source of results: 'gseapy', 'goatools', 'custom'
    sort_by : str, default 'Adjusted P-value'
        Column to sort by
    ascending : bool, default True
        Sort order
    top_n : int, optional
        Return only top N results
    significance_cutoff : float, default 0.05
        Filter by significance cutoff
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    pd.DataFrame
        Formatted results DataFrame
    """
    if isinstance(results, pd.DataFrame):
        df = results.copy()
    elif isinstance(results, list):
        # Convert list to DataFrame (assuming goatools format)
        if source == 'goatools':
            data = []
            for result in results:
                data.append({
                    'Term': getattr(result, 'name', ''),
                    'GO_ID': getattr(result, 'GO', ''),
                    'P-value': getattr(result, 'p_uncorrected', np.nan),
                    'Adjusted P-value': getattr(result, 'p_fdr_bh', np.nan),
                    'Enrichment': getattr(result, 'enrichment', np.nan),
                    'Study_count': getattr(result, 'study_count', 0),
                    'Genes': ';'.join(getattr(result, 'study_items', []))
                })
            df = pd.DataFrame(data)
        else:
            raise ValueError(f"Cannot convert list results from source: {source}")
    else:
        raise ValueError(f"Unsupported results type: {type(results)}")
    
    # Filter by significance if column exists
    adj_p_col = None
    for col in ['Adjusted P-value', 'Adjusted_P_value', 'FDR', 'q_value']:
        if col in df.columns:
            adj_p_col = col
            break
    
    if adj_p_col is not None:
        df = df[df[adj_p_col] < significance_cutoff].copy()
        if verbose:
            print(f"Filtered to {len(df)} significant results (FDR < {significance_cutoff})")
    
    # Sort results
    if sort_by in df.columns:
        df = df.sort_values(sort_by, ascending=ascending)
    
    # Return top N if specified
    if top_n is not None:
        df = df.head(top_n)
    
    if verbose:
        print(f"Returning {len(df)} formatted results")
    
    return df


def extract_genes_from_anndata(
    adata,
    layer: Optional[str] = None,
    gene_names_column: str = 'gene_name',
    use_raw: bool = False
) -> Tuple[pd.DataFrame, List[str]]:
    """
    Extract gene expression data and gene names from AnnData object
    
    Parameters
    ----------
    adata : anndata.AnnData
        AnnData object containing gene expression data
    layer : str, optional
        Layer to extract (if None, uses X)
    gene_names_column : str, default 'gene_name'
        Column in adata.var containing gene names
    use_raw : bool, default False
        Whether to use raw data
        
    Returns
    -------
    tuple
        (expression_dataframe, gene_names_list)
    """
    import anndata as ad
    
    if not isinstance(adata, ad.AnnData):
        raise ValueError("Input must be an AnnData object")
    
    # Get expression data
    if use_raw and adata.raw is not None:
        X = adata.raw.X
        var = adata.raw.var
    else:
        if layer is not None:
            X = adata.layers[layer]
        else:
            X = adata.X
        var = adata.var
    
    # Convert to DataFrame
    if hasattr(X, 'toarray'):  # sparse matrix
        X = X.toarray()
    
    df = pd.DataFrame(
        X.T,  # Transpose so genes are rows, samples are columns
        index=var.index,
        columns=adata.obs.index
    )
    
    # Get gene names
    if gene_names_column in var.columns:
        gene_names = var[gene_names_column].tolist()
    else:
        gene_names = var.index.tolist()
    
    return df, gene_names
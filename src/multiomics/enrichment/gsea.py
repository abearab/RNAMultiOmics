"""
Gene Set Enrichment Analysis (GSEA) wrapper functions using gseapy
"""

import pandas as pd
import numpy as np
from typing import Union, Optional, Dict, Any, List

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False


def run_gsea(
    data: pd.DataFrame,
    gene_sets: Union[str, List[str], Dict[str, List[str]]],
    cls: Union[str, List[str]],
    outdir: str = 'gsea_results',
    permutation_num: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
    seed: int = 123,
    verbose: bool = False,
    **kwargs
) -> Optional[object]:
    """
    Run Gene Set Enrichment Analysis using gseapy
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data with genes as rows and samples as columns
    gene_sets : str, list, or dict
        Gene sets database. Can be a string (e.g., 'KEGG_2016', 'GO_Biological_Process_2018'),
        a list of gene set names, or a dictionary with gene set names as keys and gene lists as values
    cls : str or list
        Sample phenotype labels. Either a file path or list of labels
    outdir : str, default 'gsea_results'
        Output directory for results
    permutation_num : int, default 1000
        Number of permutations
    min_size : int, default 15
        Minimum gene set size
    max_size : int, default 500
        Maximum gene set size
    seed : int, default 123
        Random seed for reproducibility
    verbose : bool, default False
        Verbose output
    **kwargs
        Additional arguments passed to gseapy.gsea
        
    Returns
    -------
    object or None
        GSEAPY results object if successful, None if gseapy not available
    """
    if not GSEAPY_AVAILABLE:
        raise ImportError(
            "gseapy is required for GSEA analysis. Install with: pip install gseapy"
        )
    
    if verbose:
        print(f"Running GSEA with {len(data)} genes and {data.shape[1]} samples")
        
    # Run GSEA
    gs_res = gp.gsea(
        data=data,
        gene_sets=gene_sets,
        cls=cls,
        outdir=outdir,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        verbose=verbose,
        **kwargs
    )
    
    if verbose:
        print(f"GSEA completed. Results saved to {outdir}")
        
    return gs_res


def run_prerank(
    rnk: Union[pd.Series, pd.DataFrame, str],
    gene_sets: Union[str, List[str], Dict[str, List[str]]],
    outdir: str = 'prerank_results',
    permutation_num: int = 1000,
    min_size: int = 15,
    max_size: int = 500,
    seed: int = 123,
    verbose: bool = False,
    **kwargs
) -> Optional[object]:
    """
    Run pre-ranked Gene Set Enrichment Analysis using gseapy
    
    Parameters
    ----------
    rnk : pd.Series, pd.DataFrame, or str
        Pre-ranked gene list. Can be a pandas Series with gene names as index and ranks as values,
        a DataFrame with two columns (gene, rank), or a file path to a .rnk file
    gene_sets : str, list, or dict
        Gene sets database. Can be a string (e.g., 'KEGG_2016', 'GO_Biological_Process_2018'),
        a list of gene set names, or a dictionary with gene set names as keys and gene lists as values
    outdir : str, default 'prerank_results'
        Output directory for results
    permutation_num : int, default 1000
        Number of permutations
    min_size : int, default 15
        Minimum gene set size
    max_size : int, default 500
        Maximum gene set size
    seed : int, default 123
        Random seed for reproducibility
    verbose : bool, default False
        Verbose output
    **kwargs
        Additional arguments passed to gseapy.prerank
        
    Returns
    -------
    object or None
        GSEAPY results object if successful, None if gseapy not available
    """
    if not GSEAPY_AVAILABLE:
        raise ImportError(
            "gseapy is required for pre-ranked GSEA analysis. Install with: pip install gseapy"
        )
    
    if verbose:
        print("Running pre-ranked GSEA")
        
    # Run pre-ranked GSEA
    pre_res = gp.prerank(
        rnk=rnk,
        gene_sets=gene_sets,
        outdir=outdir,
        permutation_num=permutation_num,
        min_size=min_size,
        max_size=max_size,
        seed=seed,
        verbose=verbose,
        **kwargs
    )
    
    if verbose:
        print(f"Pre-ranked GSEA completed. Results saved to {outdir}")
        
    return pre_res


def create_ranked_list(
    data: pd.DataFrame,
    group1: List[str],
    group2: List[str],
    method: str = 'log2fc',
    verbose: bool = False
) -> pd.Series:
    """
    Create a ranked gene list for pre-ranked GSEA from expression data
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data with genes as rows and samples as columns
    group1 : list
        Sample names for group 1
    group2 : list
        Sample names for group 2
    method : str, default 'log2fc'
        Method for ranking genes. Options: 'log2fc', 'fold_change', 't_test'
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    pd.Series
        Ranked gene list with gene names as index and ranking metric as values
    """
    # Subset data to specified groups
    group1_data = data[group1]
    group2_data = data[group2]
    
    if method == 'log2fc':
        # Calculate log2 fold change
        mean1 = group1_data.mean(axis=1)
        mean2 = group2_data.mean(axis=1)
        # Add small pseudocount to avoid log(0)
        log2fc = np.log2((mean1 + 1e-8) / (mean2 + 1e-8))
        ranked_genes = log2fc.sort_values(ascending=False)
        
    elif method == 'fold_change':
        # Calculate simple fold change
        mean1 = group1_data.mean(axis=1)
        mean2 = group2_data.mean(axis=1)
        fc = mean1 / (mean2 + 1e-8)
        ranked_genes = fc.sort_values(ascending=False)
        
    elif method == 't_test':
        from scipy import stats
        # Perform t-test for each gene
        t_stats = []
        for gene in data.index:
            t_stat, _ = stats.ttest_ind(group1_data.loc[gene], group2_data.loc[gene])
            t_stats.append(t_stat)
        
        t_series = pd.Series(t_stats, index=data.index)
        ranked_genes = t_series.sort_values(ascending=False)
        
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'log2fc', 'fold_change', 't_test'")
    
    if verbose:
        print(f"Created ranked gene list with {len(ranked_genes)} genes using {method}")
        
    return ranked_genes
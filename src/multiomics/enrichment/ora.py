"""
Over-Representation Analysis (ORA) using gseapy and statistical methods
"""

import pandas as pd
import numpy as np
from typing import Union, Optional, List, Dict, Any
from scipy.stats import fisher_exact, hypergeom

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False


def run_ora(
    gene_list: List[str],
    gene_sets: Union[str, List[str], Dict[str, List[str]]],
    background: Optional[List[str]] = None,
    organism: str = 'Human',
    outdir: str = 'ora_results',
    cutoff: float = 0.05,
    verbose: bool = False,
    **kwargs
) -> Optional[pd.DataFrame]:
    """
    Run Over-Representation Analysis using gseapy
    
    Parameters
    ----------
    gene_list : list
        List of gene symbols to test for enrichment
    gene_sets : str, list, or dict
        Gene sets database. Can be a string (e.g., 'KEGG_2016', 'GO_Biological_Process_2018'),
        a list of gene set names, or a dictionary with gene set names as keys and gene lists as values
    background : list, optional
        Background gene list. If None, uses all genes in the gene sets
    organism : str, default 'Human'
        Organism name (Human, Mouse, etc.)
    outdir : str, default 'ora_results'
        Output directory for results
    cutoff : float, default 0.05
        FDR cutoff for significance
    verbose : bool, default False
        Verbose output
    **kwargs
        Additional arguments passed to gseapy.enrichr
        
    Returns
    -------
    pd.DataFrame or None
        ORA results if successful, None if gseapy not available
    """
    if not GSEAPY_AVAILABLE:
        raise ImportError(
            "gseapy is required for ORA analysis. Install with: pip install gseapy"
        )
    
    if verbose:
        print(f"Running ORA for {len(gene_list)} genes")
        
    # Run enrichment analysis using enrichr (which performs ORA)
    enr = gp.enrichr(
        gene_list=gene_list,
        gene_sets=gene_sets,
        organism=organism,
        outdir=outdir,
        background=background,
        cutoff=cutoff,
        verbose=verbose,
        **kwargs
    )
    
    # Get results as DataFrame
    results = enr.results
    
    if verbose:
        n_significant = len(results[results['Adjusted P-value'] < cutoff])
        print(f"Found {n_significant} significantly enriched terms (FDR < {cutoff})")
        print(f"Results saved to {outdir}")
        
    return results


def run_fisher_exact_test(
    study_genes: List[str],
    gene_set: List[str],
    background_genes: List[str],
    alternative: str = 'greater'
) -> Dict[str, float]:
    """
    Run Fisher's exact test for gene set over-representation
    
    Parameters
    ----------
    study_genes : list
        List of study genes
    gene_set : list
        List of genes in the gene set of interest
    background_genes : list
        List of background genes
    alternative : str, default 'greater'
        Alternative hypothesis: 'two-sided', 'less', or 'greater'
        
    Returns
    -------
    dict
        Dictionary with test statistics and p-value
    """
    # Convert to sets for efficiency
    study_set = set(study_genes)
    geneset_set = set(gene_set)
    background_set = set(background_genes)
    
    # Calculate contingency table
    # genes in both study and gene set
    overlap = len(study_set & geneset_set)
    # genes in study but not in gene set
    study_not_geneset = len(study_set - geneset_set)
    # genes in gene set but not in study
    geneset_not_study = len((geneset_set & background_set) - study_set)
    # genes in neither study nor gene set
    neither = len(background_set - study_set - geneset_set)
    
    # Create contingency table
    # [[overlap, study_not_geneset],
    #  [geneset_not_study, neither]]
    contingency_table = [[overlap, study_not_geneset],
                        [geneset_not_study, neither]]
    
    # Run Fisher's exact test
    oddsratio, pvalue = fisher_exact(contingency_table, alternative=alternative)
    
    return {
        'overlap': overlap,
        'study_not_geneset': study_not_geneset,
        'geneset_not_study': geneset_not_study,
        'neither': neither,
        'odds_ratio': oddsratio,
        'p_value': pvalue,
        'enrichment_score': overlap / len(study_set) if len(study_set) > 0 else 0,
        'gene_set_size': len(geneset_set & background_set),
        'study_size': len(study_set)
    }


def run_hypergeometric_test(
    study_genes: List[str],
    gene_set: List[str],
    background_genes: List[str]
) -> Dict[str, float]:
    """
    Run hypergeometric test for gene set over-representation
    
    Parameters
    ----------
    study_genes : list
        List of study genes
    gene_set : list
        List of genes in the gene set of interest
    background_genes : list
        List of background genes
        
    Returns
    -------
    dict
        Dictionary with test statistics and p-value
    """
    # Convert to sets for efficiency
    study_set = set(study_genes)
    geneset_set = set(gene_set)
    background_set = set(background_genes)
    
    # Calculate parameters for hypergeometric test
    M = len(background_set)  # Population size
    n = len(geneset_set & background_set)  # Number of success states in population
    N = len(study_set)  # Number of draws
    k = len(study_set & geneset_set)  # Number of observed successes
    
    # Calculate p-value using hypergeometric distribution
    # P(X >= k) where X ~ Hypergeometric(M, n, N)
    pvalue = hypergeom.sf(k-1, M, n, N)
    
    # Calculate enrichment score
    expected = (N * n) / M if M > 0 else 0
    enrichment = k / expected if expected > 0 else float('inf')
    
    return {
        'overlap': k,
        'expected': expected,
        'enrichment_ratio': enrichment,
        'p_value': pvalue,
        'gene_set_size': n,
        'study_size': N,
        'population_size': M
    }


def run_custom_ora(
    study_genes: List[str],
    gene_sets: Dict[str, List[str]],
    background_genes: Optional[List[str]] = None,
    method: str = 'fisher',
    alpha: float = 0.05,
    multiple_testing: str = 'fdr_bh',
    verbose: bool = False
) -> pd.DataFrame:
    """
    Run custom Over-Representation Analysis using statistical tests
    
    Parameters
    ----------
    study_genes : list
        List of study genes
    gene_sets : dict
        Dictionary with gene set names as keys and gene lists as values
    background_genes : list, optional
        Background gene list. If None, uses union of all gene sets
    method : str, default 'fisher'
        Statistical test method: 'fisher' or 'hypergeometric'
    alpha : float, default 0.05
        Significance level
    multiple_testing : str, default 'fdr_bh'
        Multiple testing correction method
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    pd.DataFrame
        ORA results DataFrame
    """
    from statsmodels.stats.multitest import multipletests
    
    if background_genes is None:
        # Use union of all gene sets as background
        background_genes = list(set().union(*gene_sets.values()))
        if verbose:
            print(f"Using union of gene sets as background: {len(background_genes)} genes")
    
    results = []
    
    for gene_set_name, gene_set in gene_sets.items():
        if method == 'fisher':
            test_result = run_fisher_exact_test(study_genes, gene_set, background_genes)
        elif method == 'hypergeometric':
            test_result = run_hypergeometric_test(study_genes, gene_set, background_genes)
        else:
            raise ValueError(f"Unknown method: {method}. Choose from 'fisher' or 'hypergeometric'")
        
        results.append({
            'Term': gene_set_name,
            'Overlap': test_result['overlap'],
            'Gene_Set_Size': test_result.get('gene_set_size', len(set(gene_set) & set(background_genes))),
            'Study_Size': test_result.get('study_size', len(study_genes)),
            'P_value': test_result['p_value'],
            'Enrichment_Score': test_result.get('enrichment_ratio', test_result.get('enrichment_score', 0)),
            'Odds_Ratio': test_result.get('odds_ratio', None)
        })
    
    # Convert to DataFrame
    df = pd.DataFrame(results)
    
    # Apply multiple testing correction
    if len(df) > 0:
        reject, pvals_corrected, alpha_sidak, alpha_bonf = multipletests(
            df['P_value'], alpha=alpha, method=multiple_testing
        )
        df['Adjusted_P_value'] = pvals_corrected
        df['Significant'] = reject
    else:
        df['Adjusted_P_value'] = []
        df['Significant'] = []
    
    # Sort by adjusted p-value
    df = df.sort_values('Adjusted_P_value')
    
    if verbose:
        n_significant = len(df[df['Significant']])
        print(f"Found {n_significant} significantly enriched gene sets (FDR < {alpha})")
    
    return df


def format_ora_results(
    results: pd.DataFrame,
    top_n: Optional[int] = None,
    significance_cutoff: float = 0.05
) -> pd.DataFrame:
    """
    Format and filter ORA results
    
    Parameters
    ----------
    results : pd.DataFrame
        ORA results DataFrame
    top_n : int, optional
        Return only top N results by significance
    significance_cutoff : float, default 0.05
        Filter results by adjusted p-value cutoff
        
    Returns
    -------
    pd.DataFrame
        Formatted and filtered results
    """
    # Filter by significance
    if 'Adjusted P-value' in results.columns:
        significant_results = results[results['Adjusted P-value'] < significance_cutoff].copy()
    elif 'Adjusted_P_value' in results.columns:
        significant_results = results[results['Adjusted_P_value'] < significance_cutoff].copy()
    else:
        significant_results = results.copy()
    
    # Sort by adjusted p-value
    if 'Adjusted P-value' in significant_results.columns:
        significant_results = significant_results.sort_values('Adjusted P-value')
    elif 'Adjusted_P_value' in significant_results.columns:
        significant_results = significant_results.sort_values('Adjusted_P_value')
    
    # Return top N if specified
    if top_n is not None:
        significant_results = significant_results.head(top_n)
    
    return significant_results
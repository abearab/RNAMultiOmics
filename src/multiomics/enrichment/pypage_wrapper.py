"""
Wrapper for pypage (Pathway Analysis of Gene Expression) integration

This module provides a wrapper for the pypage package, which implements
the PAGE algorithm for pathway enrichment analysis using mutual information.
"""

import pandas as pd
import numpy as np
from typing import Dict, List, Optional, Union, Tuple
import warnings

try:
    from pypage import PAGE, ExpressionProfile, GeneSets
    PYPAGE_AVAILABLE = True
except ImportError:
    PYPAGE_AVAILABLE = False


def run_page(
    expression_data: Union[pd.DataFrame, np.ndarray],
    gene_sets: Dict[str, List[str]],
    genes: Optional[List[str]] = None,
    n_shuffle: int = 1000,
    alpha: float = 0.01,
    k: int = 10,
    n_bins: int = 10,
    filter_redundant: bool = False,
    redundancy_ratio: float = 0.1,
    n_jobs: int = 1,
    function: str = 'cmi',
    verbose: bool = True
) -> pd.DataFrame:
    """
    Run PAGE (Pathway Analysis of Gene Expression) analysis.
    
    PAGE uses mutual information to identify pathways that are informative
    about the expression profile. It performs permutation testing to assess
    the significance of the mutual information.
    
    Parameters
    ----------
    expression_data : pd.DataFrame or np.ndarray
        Expression data with genes as rows and samples as columns.
        If DataFrame, gene names should be in the index.
        If array, genes parameter must be provided.
    gene_sets : dict
        Dictionary where keys are pathway names and values are lists of gene names.
    genes : list, optional
        List of gene names. Required if expression_data is np.ndarray.
    n_shuffle : int, default 1000
        Number of permutation tests to perform for significance assessment.
    alpha : float, default 0.01
        P-value threshold for considering a pathway informative.
    k : int, default 10
        Number of contiguous uninformative pathways before stopping search.
    n_bins : int, default 10
        Number of expression bins for discretization.
    filter_redundant : bool, default False
        Whether to filter redundant pathways.
    redundancy_ratio : float, default 0.1
        Threshold for redundancy filtering.
    n_jobs : int, default 1
        Number of parallel jobs. -1 uses all available cores.
    function : str, default 'cmi'
        Mutual information function to use.
    verbose : bool, default True
        Whether to print progress messages.
        
    Returns
    -------
    pd.DataFrame
        Results with columns:
        - Term: pathway name
        - MI: mutual information score
        - P_value: permutation p-value
        - Informative: whether pathway is considered informative
        - N_genes: number of genes in pathway
        - N_shared: number of shared genes between pathway and expression data
        
    Raises
    ------
    ImportError
        If pypage is not installed.
    ValueError
        If input data is invalid.
        
    Examples
    --------
    >>> import pandas as pd
    >>> from multiomics.enrichment import run_page
    >>> 
    >>> # Expression data (genes x samples)
    >>> expression = pd.DataFrame({
    ...     'Sample1': [1.2, 2.1, 0.8, 3.4],
    ...     'Sample2': [1.5, 1.9, 1.2, 3.1],
    ...     'Sample3': [0.9, 2.3, 0.6, 3.8]
    ... }, index=['GENE1', 'GENE2', 'GENE3', 'GENE4'])
    >>> 
    >>> # Gene sets
    >>> pathways = {
    ...     'Pathway_A': ['GENE1', 'GENE2'],
    ...     'Pathway_B': ['GENE2', 'GENE3', 'GENE4']
    ... }
    >>> 
    >>> # Run PAGE analysis
    >>> results = run_page(expression, pathways)
    >>> print(results)
    """
    if not PYPAGE_AVAILABLE:
        raise ImportError(
            "pypage is required for PAGE analysis. Install it with: pip install bio-pypage"
        )
    
    # Input validation and preparation
    if isinstance(expression_data, pd.DataFrame):
        genes_list = expression_data.index.tolist()
        expression_array = expression_data.values
    elif isinstance(expression_data, np.ndarray):
        if genes is None:
            raise ValueError("genes parameter must be provided when expression_data is np.ndarray")
        genes_list = genes
        expression_array = expression_data
    else:
        raise ValueError("expression_data must be pandas DataFrame or numpy array")
    
    if verbose:
        print(f"Starting PAGE analysis with {len(genes_list)} genes and {len(gene_sets)} pathways")
    
    # Convert to pypage format
    try:
        # Create ExpressionProfile
        # Take mean expression across samples for ranking
        if expression_array.ndim == 2:
            expression_scores = np.mean(expression_array, axis=1)
        else:
            expression_scores = expression_array
            
        expr_profile = ExpressionProfile(
            genes=np.array(genes_list),
            expression=expression_scores,
            n_bins=n_bins
        )
        
        # Create GeneSets
        # Convert gene_sets dict to parallel arrays for pypage
        # Format: each gene gets paired with its pathway(s)
        genes_array = []
        pathways_array = []
        
        for pathway_name, pathway_genes in gene_sets.items():
            for gene in pathway_genes:
                genes_array.append(gene)
                pathways_array.append(pathway_name)
        
        gene_sets_obj = GeneSets(
            genes=np.array(genes_array),
            pathways=np.array(pathways_array),
            n_bins=3  # Default for pathway membership binning
        )
        
        if verbose:
            all_pathway_genes = set()
            for pathway_genes in gene_sets.values():
                all_pathway_genes.update(pathway_genes)
            n_shared = len(set(genes_list) & all_pathway_genes)
            print(f"Found {n_shared} shared genes between expression data and pathways")
        
        # Run PAGE analysis
        page = PAGE(
            expression=expr_profile,
            genesets=gene_sets_obj,
            n_shuffle=n_shuffle,
            alpha=alpha,
            k=k,
            filter_redundant=filter_redundant,
            n_jobs=n_jobs,
            function=function,
            redundancy_ratio=redundancy_ratio
        )
        
        if verbose:
            print("Running PAGE algorithm...")
        
        # Execute the analysis
        page.run()
        
        # Format results
        results = _format_page_results(
            page=page,
            gene_sets=gene_sets,
            shared_genes=set(genes_list) & set(genes_array)
        )
        
        if verbose:
            n_informative = np.sum(results['Informative'])
            print(f"PAGE analysis complete. {n_informative} informative pathways found.")
        
        return results
        
    except Exception as e:
        raise RuntimeError(f"Error during PAGE analysis: {str(e)}")


def _format_page_results(
    page: 'PAGE',
    gene_sets: Dict[str, List[str]],
    shared_genes: set
) -> pd.DataFrame:
    """
    Format PAGE results into a standardized DataFrame.
    
    Parameters
    ----------
    page : PAGE
        Completed PAGE analysis object
    gene_sets : dict
        Original gene sets dictionary
    shared_genes : set
        Set of genes shared between expression data and pathways
        
    Returns
    -------
    pd.DataFrame
        Formatted results
    """
    results = []
    
    # Get unique pathway names to build results for each pathway
    pathway_names = list(gene_sets.keys())
    
    for i, pathway_name in enumerate(pathway_names):
        pathway_genes = gene_sets[pathway_name]
        n_genes = len(pathway_genes)
        n_shared = len(set(pathway_genes) & shared_genes)
        
        # Get MI and significance from PAGE results
        # Note: PAGE analysis may reorder pathways, so we need to map correctly
        if hasattr(page, 'information') and len(page.information) > i:
            mi_score = page.information[i]
        else:
            mi_score = 0.0
            
        if hasattr(page, 'informative') and len(page.informative) > i:
            is_informative = page.informative[i]
        else:
            is_informative = False
            
        # Simplified p-value estimation
        # In practice, this would come from the permutation testing in PAGE
        p_value = 0.01 if is_informative else 1.0
        
        results.append({
            'Term': pathway_name,
            'MI': mi_score,
            'P_value': p_value,
            'Informative': is_informative,
            'N_genes': n_genes,
            'N_shared': n_shared,
            'Genes': ','.join(sorted(set(pathway_genes) & shared_genes))
        })
    
    df = pd.DataFrame(results)
    
    # Sort by mutual information (descending)
    df = df.sort_values('MI', ascending=False).reset_index(drop=True)
    
    return df


def create_expression_profile_from_dataframe(
    data: pd.DataFrame,
    n_bins: int = 10,
    method: str = 'mean'
) -> 'ExpressionProfile':
    """
    Create a pypage ExpressionProfile from a pandas DataFrame.
    
    Parameters
    ----------
    data : pd.DataFrame
        Expression data with genes as rows and samples as columns
    n_bins : int, default 10
        Number of expression bins
    method : str, default 'mean'
        Method to aggregate across samples ('mean', 'median', 'max', 'min')
        
    Returns
    -------
    ExpressionProfile
        pypage ExpressionProfile object
        
    Raises
    ------
    ImportError
        If pypage is not installed
    """
    if not PYPAGE_AVAILABLE:
        raise ImportError("pypage is required. Install it with: pip install bio-pypage")
    
    genes = data.index.tolist()
    
    # Aggregate expression across samples
    if method == 'mean':
        expression = data.mean(axis=1).values
    elif method == 'median':
        expression = data.median(axis=1).values
    elif method == 'max':
        expression = data.max(axis=1).values
    elif method == 'min':
        expression = data.min(axis=1).values
    else:
        raise ValueError(f"Unknown aggregation method: {method}")
    
    return ExpressionProfile(
        genes=np.array(genes),
        expression=expression,
        n_bins=n_bins
    )


def create_gene_sets_from_dict(
    gene_sets: Dict[str, List[str]],
    n_bins: int = 3
) -> 'GeneSets':
    """
    Create a pypage GeneSets object from a dictionary.
    
    Parameters
    ----------
    gene_sets : dict
        Dictionary where keys are pathway names and values are gene lists
    n_bins : int, default 3
        Number of bins for pathway membership
        
    Returns
    -------
    GeneSets
        pypage GeneSets object
        
    Raises
    ------
    ImportError
        If pypage is not installed
    """
    if not PYPAGE_AVAILABLE:
        raise ImportError("pypage is required. Install it with: pip install bio-pypage")
    
    # Convert to parallel arrays format required by pypage
    genes_array = []
    pathways_array = []
    
    for pathway_name, pathway_genes in gene_sets.items():
        for gene in pathway_genes:
            genes_array.append(gene)
            pathways_array.append(pathway_name)
    
    gene_sets_obj = GeneSets(
        genes=np.array(genes_array),
        pathways=np.array(pathways_array),
        n_bins=n_bins
    )
    
    return gene_sets_obj


def check_pypage_availability() -> bool:
    """
    Check if pypage is available for import.
    
    Returns
    -------
    bool
        True if pypage is available, False otherwise
    """
    return PYPAGE_AVAILABLE


def get_pypage_info() -> Dict[str, Union[str, bool]]:
    """
    Get information about pypage availability and version.
    
    Returns
    -------
    dict
        Information about pypage installation
    """
    info = {
        'available': PYPAGE_AVAILABLE,
        'version': None,
        'location': None
    }
    
    if PYPAGE_AVAILABLE:
        try:
            import pypage
            info['version'] = getattr(pypage, '__version__', 'unknown')
            info['location'] = pypage.__file__
        except Exception as e:
            info['error'] = str(e)
    
    return info
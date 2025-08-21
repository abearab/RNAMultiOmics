"""
Gene Ontology (GO) enrichment analysis using goatools and gseapy
"""

import pandas as pd
import numpy as np
from typing import Union, Optional, List, Set, Dict, Any

try:
    import gseapy as gp
    GSEAPY_AVAILABLE = True
except ImportError:
    GSEAPY_AVAILABLE = False

try:
    from goatools.base import download_go_basic_obo
    from goatools.base import download_ncbi_associations
    from goatools.obo_parser import GODag
    from goatools.anno.genetogo_reader import Gene2GoReader
    from goatools.goea.go_enrichment_ns import GOEnrichmentStudyNS
    GOATOOLS_AVAILABLE = True
except ImportError:
    GOATOOLS_AVAILABLE = False


def run_go_enrichment_gseapy(
    gene_list: List[str],
    organism: str = 'Human',
    gene_sets: Union[str, List[str]] = ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023'],
    background: Optional[List[str]] = None,
    outdir: str = 'go_enrichment_results',
    cutoff: float = 0.05,
    verbose: bool = False,
    **kwargs
) -> Optional[pd.DataFrame]:
    """
    Run GO enrichment analysis using gseapy
    
    Parameters
    ----------
    gene_list : list
        List of gene symbols to test for enrichment
    organism : str, default 'Human'
        Organism name (Human, Mouse, etc.)
    gene_sets : str or list, default ['GO_Biological_Process_2023', 'GO_Cellular_Component_2023', 'GO_Molecular_Function_2023']
        GO gene sets to use for enrichment
    background : list, optional
        Background gene list. If None, uses all genes in the gene sets
    outdir : str, default 'go_enrichment_results'
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
        Enrichment results if successful, None if gseapy not available
    """
    if not GSEAPY_AVAILABLE:
        raise ImportError(
            "gseapy is required for GO enrichment analysis. Install with: pip install gseapy"
        )
    
    if verbose:
        print(f"Running GO enrichment analysis for {len(gene_list)} genes")
        
    # Run enrichment analysis
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
        print(f"Found {n_significant} significantly enriched GO terms (FDR < {cutoff})")
        print(f"Results saved to {outdir}")
        
    return results


def run_go_enrichment_goatools(
    study_genes: Set[str],
    population_genes: Set[str],
    obo_file: Optional[str] = None,
    association_file: Optional[str] = None,
    taxid: int = 9606,  # Human
    alpha: float = 0.05,
    method: str = 'fdr_bh',
    verbose: bool = False
) -> Optional[List]:
    """
    Run GO enrichment analysis using goatools
    
    Parameters
    ----------
    study_genes : set
        Set of study genes (gene symbols or IDs)
    population_genes : set
        Set of population/background genes
    obo_file : str, optional
        Path to GO OBO file. If None, will download automatically
    association_file : str, optional
        Path to gene association file. If None, will download for specified taxid
    taxid : int, default 9606
        NCBI taxonomy ID (9606 for human, 10090 for mouse)
    alpha : float, default 0.05
        Significance level
    method : str, default 'fdr_bh'
        Multiple testing correction method
    verbose : bool, default False
        Verbose output
        
    Returns
    -------
    list or None
        GO enrichment results if successful, None if goatools not available
    """
    if not GOATOOLS_AVAILABLE:
        raise ImportError(
            "goatools is required for GO enrichment analysis. Install with: pip install goatools"
        )
    
    if verbose:
        print(f"Running GO enrichment analysis with goatools")
        print(f"Study genes: {len(study_genes)}, Population genes: {len(population_genes)}")
    
    # Download GO OBO file if not provided
    if obo_file is None:
        obo_file = download_go_basic_obo()
        if verbose:
            print(f"Downloaded GO OBO file: {obo_file}")
    
    # Load GO DAG
    obodag = GODag(obo_file)
    
    # Download gene associations if not provided
    if association_file is None:
        association_file = download_ncbi_associations()
        if verbose:
            print(f"Downloaded gene associations: {association_file}")
    
    # Load gene-to-GO associations
    objanno = Gene2GoReader(association_file, taxids=[taxid])
    ns2assoc = objanno.get_ns2assc()
    
    # Run GO enrichment analysis
    goeaobj = GOEnrichmentStudyNS(
        population_genes,
        ns2assoc,
        obodag,
        propagate_counts=False,
        alpha=alpha,
        methods=[method]
    )
    
    goea_results_all = goeaobj.run_study(study_genes)
    goea_results_sig = [r for r in goea_results_all if r.p_fdr_bh < alpha]
    
    if verbose:
        print(f"Found {len(goea_results_sig)} significantly enriched GO terms (FDR < {alpha})")
    
    return goea_results_sig


def run_go_enrichment(
    gene_list: List[str],
    method: str = 'gseapy',
    background: Optional[List[str]] = None,
    organism: str = 'Human',
    cutoff: float = 0.05,
    outdir: str = 'go_enrichment_results',
    verbose: bool = False,
    **kwargs
) -> Optional[Union[pd.DataFrame, List]]:
    """
    Run GO enrichment analysis using either gseapy or goatools
    
    Parameters
    ----------
    gene_list : list
        List of gene symbols to test for enrichment
    method : str, default 'gseapy'
        Method to use: 'gseapy' or 'goatools'
    background : list, optional
        Background gene list
    organism : str, default 'Human'
        Organism name
    cutoff : float, default 0.05
        Significance cutoff
    outdir : str, default 'go_enrichment_results'
        Output directory for results
    verbose : bool, default False
        Verbose output
    **kwargs
        Additional arguments passed to the underlying function
        
    Returns
    -------
    pd.DataFrame or list or None
        Enrichment results
    """
    if method == 'gseapy':
        return run_go_enrichment_gseapy(
            gene_list=gene_list,
            organism=organism,
            background=background,
            outdir=outdir,
            cutoff=cutoff,
            verbose=verbose,
            **kwargs
        )
    elif method == 'goatools':
        study_genes = set(gene_list)
        if background is not None:
            population_genes = set(background)
        else:
            # Use study genes as population if no background provided
            population_genes = study_genes
            if verbose:
                print("Warning: No background provided, using study genes as population")
        
        # Map organism to taxid
        organism_to_taxid = {
            'Human': 9606,
            'Mouse': 10090,
            'Rat': 10116,
            'Zebrafish': 7955,
            'Fly': 7227,
            'Worm': 6239,
            'Yeast': 559292
        }
        taxid = organism_to_taxid.get(organism, 9606)
        
        return run_go_enrichment_goatools(
            study_genes=study_genes,
            population_genes=population_genes,
            taxid=taxid,
            alpha=cutoff,
            verbose=verbose,
            **kwargs
        )
    else:
        raise ValueError(f"Unknown method: {method}. Choose from 'gseapy' or 'goatools'")


def format_go_results(
    results: Union[pd.DataFrame, List],
    method: str = 'gseapy',
    top_n: Optional[int] = None
) -> pd.DataFrame:
    """
    Format GO enrichment results into a standardized DataFrame
    
    Parameters
    ----------
    results : pd.DataFrame or list
        GO enrichment results from gseapy or goatools
    method : str, default 'gseapy'
        Method used to generate results
    top_n : int, optional
        Return only top N results by significance
        
    Returns
    -------
    pd.DataFrame
        Formatted results DataFrame
    """
    if method == 'gseapy':
        # Results are already in DataFrame format
        df = results.copy()
        if top_n is not None:
            df = df.head(top_n)
    
    elif method == 'goatools':
        # Convert goatools results to DataFrame
        data = []
        for goea_result in results:
            data.append({
                'Term': goea_result.name,
                'GO_ID': goea_result.GO,
                'P-value': goea_result.p_uncorrected,
                'Adjusted P-value': goea_result.p_fdr_bh,
                'Overlap': f"{goea_result.study_count}/{goea_result.study_items}",
                'Genes': ';'.join(goea_result.study_items),
                'Enrichment_Score': goea_result.enrichment,
                'Namespace': goea_result.NS
            })
        
        df = pd.DataFrame(data)
        if top_n is not None:
            df = df.head(top_n)
    
    else:
        raise ValueError(f"Unknown method: {method}")
    
    return df
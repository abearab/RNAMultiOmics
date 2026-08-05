"""
Gene Set Enrichment Analysis module for RNAMultiOmics

This module provides functions for pathway and gene-set enrichment analysis,
including GSEA, GO enrichment, over-representation analysis, and PAGE analysis.
"""

from .gsea import run_gsea, run_prerank
from .go_analysis import run_go_enrichment
from .ora import run_ora
from ._utils import prepare_gene_list, format_results
from .pypage_wrapper import (
    run_page, 
    create_expression_profile_from_dataframe,
    create_gene_sets_from_dict,
    check_pypage_availability,
    get_pypage_info
)

__all__ = [
    'run_gsea',
    'run_prerank', 
    'run_go_enrichment',
    'run_ora',
    'prepare_gene_list',
    'format_results',
    'run_page',
    'create_expression_profile_from_dataframe',
    'create_gene_sets_from_dict',
    'check_pypage_availability',
    'get_pypage_info'
]
import os
import polars as pl
import biobear as bb
import anndata as ad
import pandas as pd
from glob import glob

import scanpy as sc
from pydeseq2.dds import DeseqDataSet
from pydeseq2.default_inference import DefaultInference
from pydeseq2.ds import DeseqStats

from ..utils._squab_feat_cnt import load_squab_counts


def preprocess_counts(adata, method="deseq2", layer="raw_counts", min_counts=1, n_cpus=8, verbose=True):
    """Normalize counts in anndata object using specified method
    """

    adata.X = adata.layers[layer].copy()
    if min_counts != None or min_counts > 0:
        if verbose: print(f'Filtering genes with no counts / low counts (min_counts={min_counts})...')
        sc.pp.filter_genes(adata, min_counts=min_counts) # filter out genes with no counts

    if method == "deseq2":
        if verbose: print('Normalizing counts using DESeq2 method...')

        inference = DefaultInference(n_cpus=n_cpus)

        if verbose: print(f'\tcreating `dds` object...')
        
        dds = DeseqDataSet(
            counts=adata.to_df(),
            metadata=adata.obs,
            design_factors='',
            refit_cooks=True,
            inference=inference,
            quiet=True
        )
        
        dds.var = adata.var.copy()

        dds.deseq2()

        if verbose: print(f'\tadding normalized counts to anndata object...')

        adata.obs['size_factors'] = dds.obs['size_factors']
        adata.layers['deseq2_normalized'] = dds.to_df().values

    else:
        raise ValueError(f"Normalization method '{method}' not recognized. Please use 'deseq2'.")

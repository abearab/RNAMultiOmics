import os
import polars as pl
import pandas as pd
from glob import glob

import anndata as ad
import biobear as bb

from ._annotations import load_gtf, create_gene2name


def load_squab_counts(squab_dir, feature2name, feature_id_attr='gene_id', verbose=False):
    """Read squab output files and return anndata object
    load squab output files from `squab_dir` for raw counts
    
    Parameters
    ----------
    squab_dir : str
        Path to the directory containing squab output files.
    feature2name : pl.DataFrame
        A Polars DataFrame containing the mapping of feature IDs to feature names.
    feature_id_attr : str, optional
        The attribute name for the feature ID (default is 'gene_id').
    """

    # Load raw counts
    if verbose: print('Loading raw counts...')
    raw_counts = _read_squab_files(squab_dir, ".counts.tsv")
    raw_counts.index.name = feature_id_attr
    
    # Create anndata object
    if verbose: print('Creating anndata object...')
    adata = ad.AnnData(X=raw_counts.T)
    adata.layers["raw_counts"] = raw_counts.values.T

    adata.var = feature2name.set_index(feature_id_attr).loc[adata.var.index,:]

    if verbose: print(f'counts for {adata.shape[1]} features and {adata.shape[0]} samples loaded successfully.')

    return adata


def _read_squab_files(squab_dir, suffix):
    files = glob(os.path.join(squab_dir, f"*{suffix}"))
    dfs = []
    for f in files:
        sample_id = os.path.basename(f).replace(suffix, "")
        df = pd.read_csv(f, sep="\t", index_col=0, header=None)
        df = df.loc[~df.index.str.contains("__", regex=False, na=False), :]
        df.columns = [sample_id]
        dfs.append(df)

    out = pd.concat(dfs, axis=1)

    return out

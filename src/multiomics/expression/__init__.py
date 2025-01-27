import os
import polars as pl
import biobear as bb
import anndata as ad
import pandas as pd
from glob import glob
from pytximport import tximport

from .._annotations import load_gtf, create_gene2name, create_tx2gene


def load_salmon_quants(quants_dir, pattern, GTF, verbose=False):
    """Load salmon quantification files and return anndata object
    load salmon quantification files from `quants_dir` with `pattern` for quant.sf files
    """
    if type(GTF) == str:
        gtf_df = load_gtf(GTF, output_type='pl', verbose=verbose)
    elif type(GTF) == pl.DataFrame or type(GTF) == pd.DataFrame:
        gtf_df = GTF
        if verbose: print('Using provided GTF DataFrame.')
    
    tx2gene = create_tx2gene(gtf_df, verbose=verbose)
    gene2name = create_gene2name(gtf_df, verbose=verbose)
    
    rnaseq_data = tximport(
        glob(f'{quants_dir}/{pattern}/quant.sf'),
        "salmon",
        tx2gene,
    )
    
    rnaseq_data.obs.index = [x.replace(f'{quants_dir}/','').replace('/quant.sf','') for x in rnaseq_data.obs.index]
    rnaseq_data.var = gene2name.set_index('gene_id').loc[rnaseq_data.var.index,:]

    return rnaseq_data


def load_squab_counts(squab_dir, GTF, verbose=False):
    """Read squab output files and return anndata object
    load squab output files from `squab_dir` for raw counts, FPKM and TPM
    """
    if type(GTF) == str:
        gtf_df = load_gtf(GTF, output_type='pl', verbose=verbose)
    elif type(GTF) == pl.DataFrame or type(GTF) == pd.DataFrame:
        gtf_df = GTF
        if verbose: print('Using provided GTF DataFrame.')
    
    gene2name = create_gene2name(gtf_df, verbose=verbose)

    # Load raw counts
    if verbose: print('Loading raw counts...')
    raw_counts = _read_squab_files(squab_dir, ".counts.tsv", index_col=0, header=None, comment="_")
    
    # Load FPKM
    if verbose: print('Loading FPKM normalized counts...')
    fpkm = _read_squab_files(squab_dir, ".counts.fpkm.tsv", index_col=0, header=None, skiprows=3)

    # Load TPM
    if verbose: print('Loading TPM normalized counts...')
    tpm = _read_squab_files(squab_dir, ".counts.tpm.tsv", index_col=0, header=None, skiprows=3)

    # Create anndata object
    if verbose: print('Creating anndata object...')
    adata = ad.AnnData(X=raw_counts.T)
    adata.layers["fpkm"] = fpkm.values.T
    adata.layers["tpm"] = tpm.values.T

    adata.var = gene2name.set_index('gene_id').loc[adata.var.index,:]

    if verbose: print(f'counts for {adata.shape[1]} features and {adata.shape[0]} samples loaded successfully.')

    return adata


def _read_squab_files(squab_dir, suffix, **kwargs):
    files = glob(os.path.join(squab_dir, f"*{suffix}"))
    dfs = []
    for f in files:
        sample_id = os.path.basename(f).replace(suffix, "")
        df = pd.read_csv(f, sep="\t", **kwargs)
        df.columns = [sample_id]
        df.index.name = "gene_id"
        dfs.append(df)

    out = pd.concat(dfs, axis=1)

    return out

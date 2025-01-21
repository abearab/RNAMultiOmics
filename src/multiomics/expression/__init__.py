import os
import polars as pl
import biobear as bb
import anndata as ad
import pandas as pd
from glob import glob
from pytximport import tximport


def load_gtf(GTF, verbose=False):
    session = bb.new_session()
    gtf_df = session.read_gtf_file(GTF).to_polars()

    if verbose: print('GTF file loaded successfully.')
    
    return gtf_df


def create_tx2gene(gtf_df, verbose=False):
    # Filter for transcripts and create tx2gene DataFrame
    tx2gene = gtf_df.filter(
        pl.col('type') == 'transcript'
    ).with_columns([
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, 'transcript_id'), return_dtype=pl.Utf8).alias('transcript_id'),
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, 'gene_id'), return_dtype=pl.Utf8).alias('gene_id')
        # pl.col('attributes').map_elements(lambda x: extract_attribute(x, 'gene_name'), return_dtype=pl.Utf8).alias('gene_name')
    ]).select(['transcript_id','gene_id']).to_pandas()

    if verbose: print('tx2gene mapping created successfully.')

    return tx2gene


def create_gene2name(gtf_df, verbose=False):
    # Filter for genes and create gene2name DataFrame
    gene2name = gtf_df.filter(
        pl.col('type') == 'gene'
    ).with_columns([
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, 'gene_id'), return_dtype=pl.Utf8).alias('gene_id'),
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, 'gene_name'), return_dtype=pl.Utf8).alias('gene_name')
    ]).select(['gene_id','gene_name']).to_pandas()

    if verbose: print('gene2name mapping created successfully.')

    return gene2name


def load_salmon_quants(quants_dir, pattern, GTF, verbose=False):
    if type(GTF) == str:
        gtf_df = load_gtf(GTF, verbose)
    elif type(GTF) == pl.DataFrame or type(GTF) == pd.DataFrame:
        gtf_df = GTF
        if verbose: print('Using provided GTF DataFrame.')
    
    tx2gene = create_tx2gene(gtf_df, verbose)
    gene2name = create_gene2name(gtf_df, verbose)
    
    rnaseq_data = tximport(
        glob(f'{quants_dir}/{pattern}/quant.sf'),
        "salmon",
        tx2gene,
    )
    
    rnaseq_data.obs.index = [x.replace(f'{quants_dir}/','').replace('/quant.sf','') for x in rnaseq_data.obs.index]
    rnaseq_data.var = gene2name.set_index('gene_id').loc[rnaseq_data.var.index,:]

    return rnaseq_data


def load_squab_counts(squab_dir, verbose=False):
    """Read squab output files and return anndata object
    load squab output files from `squab_dir` for raw counts, FPKM and TPM
    """
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


def _extract_attribute(attribute, key):
    out = dict([(d['key'],d['value']) for d in list(attribute)])

    return out[key]
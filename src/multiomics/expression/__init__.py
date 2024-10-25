import polars as pl
import biobear as bb

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


def _extract_attribute(attribute, key):
    out = dict([(d['key'],d['value']) for d in list(attribute)])

    return out[key]
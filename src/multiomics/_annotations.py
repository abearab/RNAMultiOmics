import polars as pl
import biobear as bb
from genomicranges import GenomicRanges


def load_gtf(GTF, output_type='GenomicRanges', verbose=False):
    session = bb.new_session()
    gtf_df = session.read_gtf_file(GTF).to_polars()

    if verbose: print('GTF file loaded successfully.')

    if output_type == 'Polars' or output_type == 'pl':
        return gtf_df
    
    elif output_type == 'GenomicRanges' or output_type == 'GRanges' or output_type == 'gr':
        # rename columns: start -> starts, end -> ends, seqname -> seqnames
        gtf_gr = GenomicRanges.from_polars(
            gtf_df.rename({
                'start': 'starts', 'end': 'ends', 'seqname': 'seqnames'
            })
        )
        if verbose: print('Data converted to GenomicRanges object.')
        return gtf_gr
    
    else:
        raise ValueError('Invalid output_type. Please use "Polars"/"pl" or "GenomicRanges"/"GRanges"/"gr".')


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


def _extract_attribute(attribute, key):
    out = dict([(d['key'],d['value']) for d in list(attribute)])

    return out[key]

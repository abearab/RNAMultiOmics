import polars as pl
import biobear as bb
from genomicranges import GenomicRanges


def load_gff3(GFF, output_type='GenomicRanges', verbose=False):

    if GFF.endswith('gff3.gz') or GFF.endswith('gff3'):
        if verbose: print('Loading GFF3 file...')
    else:
        raise ValueError('Invalid GFF3 file. Please provide a valid GFF3 file with .gff3 or .gff3.gz extension.')

    # Load into Polars, skipping comments
    df = pl.read_csv(
        source=GFF,
        separator="\t",
        has_header=False,
        comment_prefix="#",
    ).rename({
        "column_1": "seqid",
        "column_2": "source",
        "column_3": "type",
        "column_4": "start",
        "column_5": "end",
        "column_6": "score",
        "column_7": "strand",
        "column_8": "phase",
        "column_9": "attributes",
    })

    # change attributes column to list of key-value pairs and set the column back in the DataFrame
    df = df.with_columns(
        pl.col("attributes")
        .str.split(";")
        .list.eval(
            pl.element()
            .str.split_exact("=", 1)
            .struct.rename_fields(["key", "value"])
        )
    )

    if verbose: print('GFF file loaded successfully.')

    if output_type == 'Polars' or output_type == 'pl':
        return df

    elif output_type == 'GenomicRanges' or output_type == 'GRanges' or output_type == 'gr':
        # rename columns: start -> starts, end -> ends, seqid -> seqnames
        gff_gr = GenomicRanges.from_polars(
            df.rename({
                'start': 'starts', 'end': 'ends', 'seqid': 'seqnames'
            })
        )
        if verbose: print('Data converted to GenomicRanges object.')
        return gff_gr

    else:
        raise ValueError('Invalid output_type. Please use "Polars"/"pl" or "GenomicRanges"/"GRanges"/"gr".')


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


def create_feature2name(df, feature_type='gene', feature_id_attr='gene_id', feature_name_attr='gene_name', verbose=False):
    """
    Create a mapping DataFrame from feature IDs to feature names based on the specified feature type and attributes.

    Parameters:
    - df (pl.DataFrame): Input Polars DataFrame containing GTF/GFF3 data.
    - feature_type (str): The type of feature to filter for (default is 'gene').
    - feature_id_attr (str): The attribute name for the feature ID (default is 'gene_id').
    - feature_name_attr (str): The attribute name for the feature name (default is 'gene_name').
    - verbose (bool): If True, print a message when the mapping is created (default is False).

    Returns:
    - feature2name (pd.DataFrame): A Pandas DataFrame with two columns: 'feature_id' and 'feature_name'.
    """
    # Filter for the specified feature type and create the mapping DataFrame
    feature2name = df.filter(
        pl.col('type') == feature_type
    ).with_columns([
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, feature_id_attr), return_dtype=pl.Utf8).alias('feature_id'),
        pl.col('attributes').map_elements(
            lambda x: _extract_attribute(x, feature_name_attr), return_dtype=pl.Utf8).alias('feature_name')
    ]).select(['feature_id', 'feature_name']).to_pandas()

    if verbose:
        print(f'{feature_type.capitalize()} to name mapping created successfully.')

    return feature2name


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

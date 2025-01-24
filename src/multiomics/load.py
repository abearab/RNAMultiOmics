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

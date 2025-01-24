import biobear as bb
from genomicranges import GenomicRanges


def load_gtf(GTF, type='GenomicRanges', verbose=False):
    session = bb.new_session()
    gtf_df = session.read_gtf_file(GTF).to_polars()

    if verbose: print('GTF file loaded successfully.')

    if type == 'Polars' or type == 'pl':
        return gtf_df
    
    elif type == 'GenomicRanges' or type == 'GRanges' or type == 'gr':
        return bb.polars_to_gr(gtf_df)
    
    else:
        raise ValueError('Invalid type. Please use "Polars"/"pl" or "GenomicRanges"/"GRanges"/"gr".')

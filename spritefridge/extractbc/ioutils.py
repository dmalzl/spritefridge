import gzip
import regex
import logging

import pandas as pd

from io import BytesIO


READCONTENTS = ['name', 'seq', 'spacer', 'quals']


def initialize_stats(nbcs):
    stats = {
        'valid': 0, 
        'filtered': 0,
        **{i: 0 for i in range(nbcs + 1)}
    }
    return stats


def count_valid_bcs(bcs):
    # python sourcery for quick summing
    # implicitly converts string to bool
    return sum(bool(bc) for bc in bcs)


def open_fastq(filepath):
    handle = (
        gzip.open(filepath, 'rt') 
        if filepath.endswith('gz') 
        else open(filepath, 'r')
    )
    return handle


def get_read(fastq):
    read = {
        k: fastq.readline().rstrip()
        for k
        in READCONTENTS
    }
    
    if not all(read.values()):
        return {}

    # this is necessary to ensure the aligner does not strip the barcodes later
    read['name'] = read['name'].split(maxsplit = 1)[0]
    return read


def read_fastqs(fastq1_path, fastq2_path):
    fastq1 = open_fastq(fastq1_path)
    fastq2 = open_fastq(fastq2_path)
    read1 = get_read(fastq1)
    read2 = get_read(fastq2)
    while read1 and read2:
        yield read1, read2

        read1 = get_read(fastq1)
        read2 = get_read(fastq2)

    fastq1.close()
    fastq2.close()


def compress_read(read):
    string = '\n'.join([read[k] for k in READCONTENTS]) + '\n'
    return gzip.compress(string.encode('utf-8'))


def reads_to_byteblocks(reads):
    stats = {}
    bytestreams = dict(
        r1 = BytesIO(),
        filtered_r1 = BytesIO(),
        r2 = BytesIO(),
        filtered_r2 = BytesIO()
    )
    
    for read1, read2, bcs in reads:
        # initializing here to avoid passing number fo barcodes
        if not stats:
            nbcs = len(bcs)
            stats = initialize_stats(nbcs)

        bcs_string = '|'.join(bcs)
        read1['name'] = read1['name'] + '[' + bcs_string
        read2['name'] = read2['name'] + '[' + bcs_string
        if not all(bcs):
            stats['filtered'] += 1
            stats[count_valid_bcs(bcs)] += 1
            bytestreams['filtered_r1'].write(
                compress_read(read1)
            )
            bytestreams['filtered_r2'].write(
                compress_read(read2)
            )
            continue
        
        stats['valid'] += 1
        stats[nbcs] += 1
        bytestreams['r1'].write(
            compress_read(read1)
        )
        bytestreams['r2'].write(
            compress_read(read2)
        )

    return {k: stream.getvalue() for k, stream in bytestreams.items()}, stats


def write_byteblocks(byteblocks, outfilepaths):
    for k, block in byteblocks.items():
        outfile = outfilepaths[k]
        if not outfile:
            continue

        # need to use simple file here otherwise double compression
        with open(outfile, 'ab') as out:
            out.write(block)


def write_fastq(reads, outfilepaths):
    byteblocks, blockstats = reads_to_byteblocks(reads)
    write_byteblocks(
        byteblocks,
        outfilepaths
    )
    return blockstats


def read_barcodes(barcodes_path, allowed_mismatches):
    barcodes = pd.read_csv(
        barcodes_path,
        sep = '\t',
        header = None,
        names = ['category', 'bcname', 'bcseq']
    )
    bc_dict = {}
    max_bc_lengths = {}
    for cat, cat_barcodes in barcodes.groupby('category'):
        bc_lens = cat_barcodes.bcseq.str.len()
        max_bc_lengths[cat] = (bc_lens.min(), bc_lens.max())
        mismatches = allowed_mismatches[cat]
        bc_dict[cat] = {
            bc.bcseq: {
                'regex': regex.compile('(' + bc.bcseq + '){s<=' + str(mismatches) + '}'), 
                'name': bc.bcname
            }
            for _, bc
            in cat_barcodes.iterrows()
        }
    
    return bc_dict, max_bc_lengths


def write_stats(stats_dict, statsfile):
    keys = stats_dict.keys()
    with open(statsfile, 'w') as file:
        file.write('\t'.join([f'nbcs_{k}' for k in keys]) + '\n')
        file.write('\t'.join(str(stats_dict[k]) for k in keys) + '\n')


def write_parallel(
    outfilepaths,
    input_queue,
    lock,
    nextractors
):
    logging.info('starting writer process')
    stats = {}
    reads_processed = 0
    while True:
        reads = input_queue.get()
        if not reads:
            nextractors -= 1

        if not nextractors:
            break
        
        reads_processed += len(reads)
        byteblocks, blockstats = reads_to_byteblocks(reads)

        if not stats:
            stats = {k: 0 for k in blockstats.keys()}

        for k, v in blockstats.items():
            stats[k] += v

        # currently this is not necessary since we only use one writerthread
        # but we leave here for the future probably
        with lock:
            write_byteblocks(
                byteblocks,
                outfilepaths
            )
        
        if reads and not reads_processed % 1e5:
            logging.info(f'processed {reads_processed} reads')

    logging.info('all reads processed, shutting down writer')
    write_stats(stats, outfilepaths['stats'])


def initialize_output(outfilepaths):
    for path in outfilepaths.values():
        if path:
            open(path, 'w')


def parse_layout(layout_string, minmax_bc_len):
    layout = []
    for bc_cat in layout_string.split('|'):
        min_len, max_len = minmax_bc_len[bc_cat]
        layout.append(
            [
                bc_cat,
                min_len,
                max_len
            ]
        )
    
    return layout


def parse_mismatches(mismatch_string):
    mismatch_dict = {}
    for cat_mismatch in mismatch_string.split(','):
        cat, mismatch = cat_mismatch.split(':')
        mismatch_dict[cat] = int(mismatch)

    return mismatch_dict


def add_spacer_info(minmax_bc_len, spacerlen):
    minmax_bc_len['SPACER'] = (spacerlen, spacerlen)
    
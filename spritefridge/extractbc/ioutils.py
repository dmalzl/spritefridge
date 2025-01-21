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
    stats['poswise'] = [0] * nbcs
    return stats


def open_fastq(filepath):
    handle = (
        gzip.open(filepath, 'rb') 
        if filepath.endswith('gz') 
        else open(filepath, 'rb')
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

    # remove any unnecessary strings from spacer
    read['spacer'] = read['spacer'][:1]

    # this is necessary to ensure the aligner does not strip the barcodes later
    read['name'] = read['name'].split(maxsplit = 1)[0]
    return read


def read_fastqs(fastq1_path, fastq2_path):
    with (
        open_fastq(fastq1_path) as fastq1,
        open_fastq(fastq2_path) as fastq2
    ):
        read1 = get_read(fastq1)
        read2 = get_read(fastq2)
        while read1 and read2:
            yield read1, read2

            read1 = get_read(fastq1)
            read2 = get_read(fastq2)


def compress_read(read):
    string = b'\n'.join([read[k] for k in READCONTENTS]) + b'\n'
    return gzip.compress(string)


def increment_poswise_counter(read_bcs, stats):
    n_valid = 0
    for i, bc in enumerate(read_bcs):
        if not bc:
            continue
        
        n_valid += 1
        stats['poswise'][i] += 1

    return n_valid
    

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

        n_valid_bcs = increment_poswise_counter(bcs, stats)
        bcs_string = b'|'.join(bcs)
        read1['name'] = read1['name'] + b'[' + bcs_string
        read2['name'] = read2['name'] + b'[' + bcs_string
        if not all(bcs):
            stats['filtered'] += 1
            stats[n_valid_bcs] += 1
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
            bytes(bc.bcseq, 'utf-8'): {
                'regex': regex.compile(
                    b'(' + bytes(bc.bcseq, 'utf-8') + 
                    b'){s<=' + 
                    bytes(str(mismatches), 'utf-8') + 
                    b'}'
                ), 
                'name': bytes(bc.bcname, 'utf-8')
            }
            for _, bc
            in cat_barcodes.iterrows()
        }
    
    return bc_dict, max_bc_lengths


def write_overall_stats(stats_dict, statsfile):
    with open(statsfile, 'w') as file:
        for i, count in stats_dict.items():
            if i == 'poswise':
                continue

            file.write(f'{i}_barcodes\t{count}' + '\n')


def write_poswise_stats(poswise_counts, bc_cats, n_reads, statsfile):
    with open(statsfile, 'w') as file:
        header = '\t'.join(
            [
                f'{bc_cat}_{i}' 
                for bc_cat, i 
                in enumerate(bc_cats)
            ]
        )
        file.write(header + '\n')
        data = '\t'.join(
            [
                '{:.2f}'.format(c/n_reads * 100) 
                for c 
                in poswise_counts
            ]
        )
        file.write(data + '\n')


def sum_stats(stats, blockstats):
    if not blockstats:
        return
    
    for k, v in blockstats.items():
        if k == 'poswise':
            continue
        
        stats[k] += v

    for i, count in enumerate(blockstats['poswise']):
        stats['poswise'][i] += count


def initialize_stats_from_blockstats(blockstats):
    stats = {}
    for k, v in blockstats.items():
        if k == 'poswise':
            stats[k] = [0] * len(v)
            continue

        stats[k] = 0

    return stats


def write_parallel(
    outfilepaths,
    input_queue,
    lock,
    nextractors,
    bc_cats
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
            stats = initialize_stats_from_blockstats(blockstats)

        sum_stats(stats, blockstats)

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
    n_reads = stats['valid'] + stats['filtered']
    write_overall_stats(stats, outfilepaths['overall_stats'])
    write_poswise_stats(
        stats['poswise'], 
        bc_cats, 
        n_reads, 
        outfilepaths['poswise_stats']
    )


def initialize_output(outfilepaths):
    for path in outfilepaths.values():
        if path:
            open(path, 'w')


def parse_layout(layout_string, minmax_bc_len, allowed_mismatches):
    layout = []
    for bc_cat in layout_string.split('|'):
        min_len, max_len = minmax_bc_len[bc_cat]
        layout.append(
            [
                bc_cat,
                min_len,
                max_len,
                # this will return None in case of unknown categories
                # in case of e.g. SPACER this is okay since we anyway skip it
                allowed_mismatches.get(bc_cat) 
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
    
from .extract import (
    extract_parallel,
    extract_barcodes
)

from .ioutils import (
    read_fastqs, 
    write_fastq,
    write_parallel,
    write_overall_stats,
    write_poswise_stats,
    initialize_stats,
    initialize_output,
    sum_stats
)

import multiprocessing as mp

import logging

BUFFERSIZE = 10000
                

def process_sequential(
    r1file, 
    r2file, 
    outfilepaths,
    bc_dicts,
    layout_r1,
    layout_r2,
    laxity
):
    stats = {}
    read_buffer = []
    reads_processed = 0
    for read1, read2 in read_fastqs(r1file, r2file):
        bcs = (
            extract_barcodes(read1, bc_dicts, layout_r1, laxity) +
            extract_barcodes(read2, bc_dicts, layout_r2, laxity)
        )

        if not stats:
            stats = initialize_stats(len(bcs))

        bcs_string = b'|'.join(bcs)
        read1['name'] = read1['name'] + bcs_string
        read2['name'] = read2['name'] + bcs_string
        read_buffer.append(
            [read1, read2, bcs]
        )

        if len(read_buffer) == BUFFERSIZE:
            blockstats = write_fastq(
                read_buffer,
                outfilepaths
            )
            sum_stats(stats, blockstats)

            read_buffer = []

        reads_processed += 1
        if not reads_processed % 1e5:
            logging.info(f'processed {reads_processed} reads')
    
    # write last reads
    blockstats = write_fastq(
        read_buffer,
        outfilepaths
    )
    sum_stats(stats, blockstats)

    bc_cats = [c[0] for c in layout_r1 + layout_r2 if c[0] != 'SPACER']
    n_reads = stats['valid'] + stats['filtered']
    write_overall_stats(stats, outfilepaths['overall_stats'])
    write_poswise_stats(
        stats['poswise'], 
        bc_cats, 
        n_reads, 
        outfilepaths['poswise_stats']
    )


def process_parallel(
    r1file, 
    r2file, 
    outfilepaths,
    bc_dicts,
    layout_r1,
    layout_r2,
    laxity,
    nprocesses
):
    # set queue size to avoid loading the full file in memory
    # queue blocks when full until space is freed
    extract_queue = mp.Queue(maxsize = 100)
    write_queue = mp.Queue(maxsize = 100)
    lock = mp.Lock()
    initialize_output(outfilepaths)

    extractor_args = (
        bc_dicts,
        layout_r1,
        layout_r2,
        laxity,
        extract_queue,
        write_queue
    )
    extractors = [
        mp.Process(
            target = extract_parallel,
            args = extractor_args
        )
        for _ in range(max(nprocesses - 2, 1))
    ]
    bc_cats = [c[0] for c in layout_r1 + layout_r2 if c[0] != 'SPACER']
    writer = mp.Process(
        target = write_parallel,
        args = (
            outfilepaths,
            write_queue,
            lock,
            len(extractors),
            bc_cats
        )
    )

    parallel_processes = [*extractors, writer]
    for process in parallel_processes:
        process.start()

    read_buffer = []
    for read1, read2 in read_fastqs(r1file, r2file):
        read_buffer.append([read1, read2])

        if len(read_buffer) == BUFFERSIZE:
            extract_queue.put(read_buffer)
            read_buffer = []

    # send last reads
    extract_queue.put(read_buffer)

    for _ in extractors:
        extract_queue.put([])

    logging.info('waiting for processes to finish')
    for process in parallel_processes:
        process.join()

    logging.info('all done!')

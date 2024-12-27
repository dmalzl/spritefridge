import logging

from .ioutils import (
    read_bam,
    clusters_to_pairs,
    write_stats
)
from collections import defaultdict


def main(args):
    clusters = defaultdict(set)
    logging.info('reading alignments and constructing clusters')
    alignments_processed = 0
    for bcs, pos in read_bam(args.bam, args.separator):
        clusters[bcs].add(pos)
        alignments_processed += 1

        if not alignments_processed % 1e5:
            logging.info(f'processed {alignments_processed} alignments')

    logging.info('finished cluster construction. writing pairs for all found sizes')
    stats = clusters_to_pairs(
        clusters, 
        args.outprefix,
        args.clustersizelow,
        args.clustersizehigh
    )
    logging.info('finished writing pairs. writing stats')
    write_stats(stats, args.outprefix + '.stats.tsv')

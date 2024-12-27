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
    alignments_processed, duplicate_alignments = 0, 0
    for bcs, pos in read_bam(args.bam, args.separator):
        if pos in clusters[bcs]:
            duplicate_alignments += 1
        
        clusters[bcs].add(pos)
        alignments_processed += 1

        if not alignments_processed % 1e5:
            logging.info(f'processed {alignments_processed} alignments')

    logging.info('finished cluster construction. write duplication stats')
    alignment_stats = {
        'unique': alignments_processed - duplicate_alignments,
        'duplicated': duplicate_alignments
    }
    write_stats(
        alignment_stats, 
        args.outprefix + '.duplicatestats.tsv'
    )
    logging.info('now writing pairs for all found sizes')
    size_stats = clusters_to_pairs(
        clusters, 
        args.outprefix,
        args.clustersizelow,
        args.clustersizehigh
    )
    logging.info('finished writing pairs. writing stats')
    write_stats(
        size_stats, 
        args.outprefix + '.sizestats.tsv'
    )

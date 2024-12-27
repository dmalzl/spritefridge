from .ioutils import (
    read_barcodes,
    parse_layout,
    parse_mismatches,
    add_spacer_info
)
from .processors import (
    process_parallel, 
    process_sequential
)
import logging


def generate_fileprefix(filename):
    prefix_tokens = filename.split('.')
    slice_end = -2 if filename.endswith('gz') else -1
    return '.'.join(prefix_tokens[:slice_end])
    

def main(args):
    logging.info('reading barcodes, layouts and setting up output files')
    allowed_mismatches = parse_mismatches(args.mismatches)
    bc_dicts, minmax_bc_len = read_barcodes(args.barcodes, allowed_mismatches)
    # add spacer info so we don't have to explicitly care about it
    add_spacer_info(minmax_bc_len, args.spacerlen)

    layout_r1 = parse_layout(args.layout1, minmax_bc_len)
    layout_r2 = parse_layout(args.layout2, minmax_bc_len)

    fileprefix = generate_fileprefix(args.output)
    outfilepaths = dict(
        r1 = args.output,
        filtered_r1 = fileprefix + '.filtered.fq.gz' if args.writefiltered else None,
        r2 = fileprefix + '.r2.fq.gz' if args.writer2 else None,
        filtered_r2 = fileprefix + '.filtered_r2.fq.gz' if args.writer2 and args.writefiltered else None,
        stats = fileprefix + '.stats.tsv'
    )

    if args.processes > 1:
        logging.info('processing in parallel')
        process_parallel(
            args.read1,
            args.read2,
            outfilepaths,
            bc_dicts,
            layout_r1,
            layout_r2,
            args.laxity,
            args.processes
        )

    else:
        logging.info('processing sequentially')
        process_sequential(
            args.read1,
            args.read2,
            outfilepaths,
            bc_dicts,
            layout_r1,
            layout_r2,
            args.laxity
        )
        
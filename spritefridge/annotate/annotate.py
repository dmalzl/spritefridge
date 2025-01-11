import logging

import pandas as pd

from cooler import Cooler, fileops
from .core import annotate_bins


def annotate_cool(coolpath, bedpaths, outprefix):
    cooler = Cooler(coolpath)

    annotated_bins = cooler.bins()[:].loc[:, ['chrom', 'start', 'end']]
    for bedpath in bedpaths:
        logging.info(f'annotating bins of {coolpath} with clusters from {bedpath}')
        tmp = annotate_bins(cooler, bedpath)

        annotated_bins = annotated_bins.merge(
            tmp,
            on = ['chrom', 'start', 'end'],
            how = 'left'
        )

    outfile = outprefix + '.tsv.gz'
    logging.info(f'writing annotated data to {outfile}')
    annotated_bins.to_csv(
        outfile,
        sep = '\t',
        index = False,
        compression = 'gzip'
    )


def annotate_mcool(mcoolpath, bedpaths, outprefix):
    for coolpath in fileops.list_coolers(mcoolpath):
        uri = mcoolpath + '::' + coolpath
        annotate_cool(
            uri,
            bedpaths,
            outprefix + coolpath.replace('/', '_')
        )


def main(args):
    if fileops.is_multires_file(args.input):
        logging.info('annotating multires cooler')
        annotate_mcool(args.input, args.bed, args.outprefix)

    else:
        logging.info('annotating single cooler')
        annotate_cool(args.input, args.bed, args.outprefix)

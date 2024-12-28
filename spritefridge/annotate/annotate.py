import logging

import pandas as pd

from cooler import Cooler, fileops
from .core import annotate_bins
from .ioutils import create_annotated_cooler


def annotate_cool(coolpath, bedpaths, outfile, mcoolfile = False):
    cooler = Cooler(coolpath)

    annotated_bins = pd.DataFrame()
    for bedpath in bedpaths:
        logging.info(f'annotating bins of {coolpath} with clusters from {bedpath}')
        tmp = annotate_bins(cooler, bedpath)
        if annotate_bins.empty:
            annotate_bins = tmp
            continue

        annotated_bins = annotated_bins.merge(
            tmp,
            on = ['chrom', 'start', 'end'],
            how = 'left'
        )

    logging.info(f'writing annotated data to {outfile}')
    create_annotated_cooler(
        coolpath,
        outfile,
        annotated_bins,
        cooler.chromnames,
        mcoolfile = mcoolfile
    )


def annotate_mcool(mcoolpath, bedpaths, outfile):
    for coolpath in fileops.list_coolers(mcoolpath):
        uri = mcoolpath + '::' + coolpath
        outuri = outfile + '::' + coolpath
        annotate_cool(
            uri,
            bedpaths,
            outuri,
            mcoolfile = True
        )


def main(args):
    for coolpath in args.input:
        if fileops.is_multires_file(coolpath):
            logging.info('annotating multires cooler')
            outfile = coolpath.replace('mcool', 'annotated.mcool')
            annotate_mcool(coolpath, args.bed, outfile)

        else:
            logging.info('annotating single cooler')
            outfile = coolpath.replace('.cool', '.annotated.cool')
            annotate_cool(coolpath, args.bed, outfile)

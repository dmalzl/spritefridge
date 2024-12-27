import logging
from .core import annotate_bins
from .ioutils import create_annotated_cooler

from cooler import Cooler


def main(args):
    bedpath = args.bed
    for coolerpath in args.input:
        cooler = Cooler(coolerpath)

        logging.info(f'annotating bins of {coolerpath} with clusters from {bedpath}')
        annotated_bins = annotate_bins(cooler, bedpath)

        outfile = coolerpath.replace('.cool', '.annotated.cool')
        logging.info(f'writing annotated data to {outfile}')
        create_annotated_cooler(
            coolerpath,
            outfile,
            annotated_bins,
            cooler.chromnames
        )

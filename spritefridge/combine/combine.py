import logging

import numpy as np

from .ioutils import read_coolers
from .core import chunked_pixels_merge
from cooler import create_cooler


def main(args):
    logging.info('reading coolers')
    coolers = read_coolers(args.input)
    logging.info('merging coolers')
    if not args.floatcounts:
        columns = ['count', 'float_count']
        dtypes = None
        logging.info('storing float counts in extra column')
    
    else:
        columns = None
        dtypes = {'count': float}
        logging.info('storing count as float')

    key = np.random.choice(list(coolers.keys()))
    bins = coolers[key].bins()[:]

    create_cooler(
        args.outfile,
        bins,
        chunked_pixels_merge(
            coolers, 
            args.nchunks, 
            not args.floatcounts
        ),
        columns = columns,
        dtypes = dtypes,
        ordered = False
    )
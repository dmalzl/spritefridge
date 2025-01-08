import logging

import numpy as np

from .ioutils import read_coolers, check_file_limit
from .core import SpriteCoolerMerger
from cooler.create import create


def main(args):
    logging.info('reading coolers')
    coolers = read_coolers(args.input)
    check_file_limit(len(coolers))
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
    assembly = coolers[key].info.get("genome-assembly", None)
    iterator = SpriteCoolerMerger(
        coolers, 
        mergebuf = args.chunksize
    )

    create(
        args.outfile,
        bins,
        iterator,
        columns=columns,
        dtypes=dtypes,
        assembly=assembly
    )

import logging
import gc

from .ioutils import (
    get_resolutons,
    check_weight,
    store_weights
)
from .core import (
    balance_ic,
    balance_kr
)
from cooler import fileops

def main(args):
    for resolution in get_resolutons(args.mcool):
        cooleruri = args.mcool + '::/resolutions/' + resolution
        outuri = args.output + '::/resolutions/' + resolution
        fileops.cp(cooleruri, outuri)

        for weight_name, per_chrom, balancefunc in zip(
            ['KR', 'perChromKR', 'ICE', 'perChromIC'],
            [False, True, False, True],
            [balance_kr, balance_kr, balance_ic, balance_ic]
        ):
            if not check_weight(cooleruri, weight_name) or args.overwrite:
                logging.info(
                    'computing {} for {}::/resolution/{}'.format(
                        weight_name, 
                        args.mcool, 
                        resolution
                    )
                )
                weights, stats = balancefunc(
                    outuri, 
                    per_chrom,
                    nproc = args.processors,
                    maxiter = args.maxiter
                )
                store_weights(
                    outuri, 
                    weights, 
                    weight_name,
                    stats,
                    overwrite = args.overwrite
                )
                del weights

            else:
                logging.info(
                    '{} weights for {}::resolution/{} already exist. Skipping!'.format(
                        weight_name, 
                        args.mcool, 
                        resolution
                    )
                )

        gc.collect()

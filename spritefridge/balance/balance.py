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

def main(args):
    for resolution in get_resolutons(args.mcool):
        cooleruri = args.mcool + '::resolutions/' + resolution

        for weight_name, per_chrom, balancetype in zip(
            ['KR', 'perChromKR'],
            [False, True],
            ['genomewide', 'per chromosome']
        ):
            if not check_weight(cooleruri, weight_name) or args.overwrite:
                logging.info(
                    'applying {} KR to {}::resolution/{}'.format(
                        balancetype, 
                        args.mcool, 
                        resolution
                    )
                )
                krweights = balance_kr(
                    cooleruri, 
                    per_chrom
                )
                store_weights(
                    cooleruri, 
                    krweights, 
                    weight_name,
                    overwrite = args.overwrite
                )
                del krweights

            else:
                logging.info(
                    '{} KR weights for {}::resolution/{} already exist. Skipping!'.format(
                        balancetype, 
                        args.mcool, 
                        resolution
                    )
                )

        for weight_name, per_chrom, balancetype in zip(
            ['ICE', 'perChromIC'], 
            [False, True], 
            ['genomewide', 'per chromosome']
        ):
            if not check_weight(cooleruri, weight_name) or args.overwrite:
                logging.info(
                    'applying {} IC to {}::resolution/{}'.format(
                        balancetype, 
                        args.mcool, 
                        resolution
                    )
                )
                icweights, stats = balance_ic(
                    cooleruri, 
                    args.processors, 
                    per_chrom
                )
                store_weights(
                    cooleruri, 
                    icweights, 
                    weight_name,
                    stats,
                    overwrite = args.overwrite
                )
                del icweights

            else:
                logging.info(
                    '{} IC weights for {}::resolution/{} already exist. Skipping!'.format(
                        balancetype, 
                        args.mcool, 
                        resolution
                    )
                )

        gc.collect()

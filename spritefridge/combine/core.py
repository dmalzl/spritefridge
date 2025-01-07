import logging

import pandas as pd
import numpy as np

from functools import partial
from cooler._reduce import CoolerMerger, merge_breakpoints


def load_and_correct(cool, size, lo, hi):
    selector = cool.pixels()
    pixels = selector[lo: hi]
    pixels['count'] *= (2/size)
    return pixels


class SpriteCoolerMerger(CoolerMerger):
    def __init__(self, coolers, mergebuf):
        super().__init__(list(coolers.values()), mergebuf)
        self.coolers = coolers

    def __iter__(self):
        # Load bin1_offset indexes lazily.
        indexes = [c.open("r")["indexes/bin1_offset"] for c in self.coolers.values()]

        # Calculate the common partition of bin1 offsets that define the epochs
        # of merging data.
        bin1_partition, cum_nrecords = merge_breakpoints(indexes, self.mergebuf)
        nrecords_per_epoch = np.diff(cum_nrecords)
        nnzs = [len(c.pixels()) for c in self.coolers.values()]

        starts = [0] * len(self.coolers)
        for i, bin1_id in enumerate(bin1_partition[1:], 2):
            stops = [index[bin1_id] for index in indexes]
            # extract, concat
            combined = pd.concat(
                [
                    load_and_correct(c, size, start, stop)
                    for (size, c), start, stop in zip(self.coolers.items(), starts, stops)
                    if (stop - start) > 0
                ],
                axis=0,
                ignore_index=True,
            )

            # sort and aggregate
            df = combined \
                .groupby(["bin1_id", "bin2_id"], sort=True) \
                .sum() \
                .reset_index()

            percent_merged = np.floor(i / len(bin1_partition) * 100)
            if not percent_merged % 10:
                logging.info(f'merged approximately {percent_merged} %')

            yield {k: v.values for k, v in df.items()}
            starts = stops

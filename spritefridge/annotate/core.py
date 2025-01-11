from pybedtools import BedTool
from .ioutils import read_bed_by_chrom

import pandas as pd

import os
import logging


def annotate_bins(cool, clusterbedfile):
    cool_chroms = set(cool.chromnames)
    annotated_bins = []
    for chrom_bed in read_bed_by_chrom(clusterbedfile):
        chrom = chrom_bed.iloc[0, 0]
        if not chrom in cool_chroms:
            continue
        
        lo, hi = cool.extent(chrom)
        chrom_bins = cool.bins()[lo:hi].loc[:, ['chrom', 'start', 'end']]
        chrom_bins['name'] = 'bin'
        a = BedTool.from_dataframe(chrom_bed)
        b = BedTool.from_dataframe(chrom_bins)
        intersect = b.intersect(a, wao = True, sorted = True).to_dataframe()
        grouped = intersect[['chrom', 'start', 'end', 'thickEnd']].groupby(['chrom', 'start', 'end'])
        annotated_bins.append(
            grouped.agg({'thickEnd': ','.join}).reset_index()
        )
    
    annotation = pd.concat(annotated_bins, ignore_index = True)
    colname = os.path.basename(clusterbedfile)
    annotation.rename(
        columns = {'thickEnd': colname},
        inplace = True
    )
    return annotation

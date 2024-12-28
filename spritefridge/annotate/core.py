from pybedtools import BedTool
from .ioutils import read_bed_by_chrom

import pandas as pd

import os


def sanitize_dtypes(annotated_bins):
    annotated_bins['chrom'] = annotated_bins.chrom.astype('category')
    annotated_bins['clusters'] = annotated_bins.clusters.astype('bytes')


def annotate_bins(cool, clusterbedfile):
    cool_chroms = set(cool.chromnames)
    annotated_bins = []
    for chrom_bed in read_bed_by_chrom(clusterbedfile):
        chrom = chrom_bed.iloc[0, 0]
        if not chrom in cool_chroms:
            continue
        
        lo, hi = cool.extent(chrom)
        chrom_bins = cool.bins()[lo:hi]
        chrom_bins['name'] = 'bin'
        a = BedTool.from_dataframe(chrom_bed)
        b = BedTool.from_dataframe(chrom_bins)
        intersect = b.intersect(a, wao = True, sorted = True).to_dataframe()
        grouped = intersect[['chrom', 'start', 'end', 'thickEnd']].groupby(['chrom', 'start', 'end'])
        annotated_bins.append(
            grouped.agg({'thickEnd': ','.join}).reset_index()
        )
    
    tmp_annotation = pd.concat(annotated_bins)
    colname = os.path.basename(clusterbedfile)
    tmp_annotation.rename(
        columns = {'thickEnd': colname},
        inplace = True
    )
    annotated_bins = cool.bins()[:].merge(
        tmp_annotation, 
        on = ['chrom', 'start', 'end'], 
        how = 'left'
    )
    annotated_bins.fillna({colname: 'None'}, inplace = True)
    sanitize_dtypes(annotated_bins)
    return annotated_bins

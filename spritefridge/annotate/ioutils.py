from cooler import fileops, create

import h5py

import pandas as pd


def to_right_dtype(x):
    try:
        return int(x)
    
    except ValueError:
        return x


def read_bed_by_chrom(bedfile):
    # this relies on the bed being sorted
    linebuffer = []
    with open(bedfile, 'r') as bed:
        current_chrom = None
        for line in bed:
            line = [
                to_right_dtype(field) 
                for field 
                in line.rstrip().split('\t')
            ]

            if not current_chrom:
                current_chrom = line[0]
            
            if current_chrom != line[0]:
                chrom_data = pd.DataFrame(
                    linebuffer,
                    columns = ['chrom', 'start', 'end', 'name']
                )
                current_chrom = line[0]
                linebuffer = [line]
                yield chrom_data
                continue

            linebuffer.append(line)


def create_annotated_cooler(source, dest, bins, chromnames, h5opts = None):
    fileops.cp(source, dest)
    h5 = h5py.File(dest, 'a')
    del h5['/']['bins']
    grp = h5['/'].create_group('bins')
    h5opts = create._create._set_h5opts(h5opts)
    create._create.write_bins(
        grp,
        bins,
        chromnames,
        h5opts
    )
    
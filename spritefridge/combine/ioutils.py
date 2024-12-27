import glob
import re
import os

import pandas as pd

from cooler import Cooler


cs_regex = re.compile('_(?P<cs>[0-9]+)_')


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


def clustersize_from_filename(filename):
    m = cs_regex.search(filename)
    return int(m.group('cs'))


def read_coolers(directory):
    coolers = {}
    for coolfile in glob.glob(directory + '/*'):
        clustersize = clustersize_from_filename(
            os.path.basename(coolfile)
        )
        coolers[clustersize] = Cooler(coolfile)
    
    return coolers

from cooler import fileops, create

import h5py
import logging

import pandas as pd


def attrs_to_dict(h5file):
    attrs = {
        key: h5file.attrs[key]
        for key
        in h5file.attrs.keys()
    }
    return attrs


def copy_attrs(source, dest):
    with (
        h5py.File(source, 'r+') as infile,
        h5py.File(dest, 'r+') as ofile
    ):
        attrs = attrs_to_dict(infile)
        ofile.attrs.update(attrs)


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
                logging.info(f'processed {current_chrom}')
                chrom_data = pd.DataFrame(
                    linebuffer,
                    columns = ['chrom', 'start', 'end', 'name']
                )
                current_chrom = line[0]
                linebuffer = [line]
                yield chrom_data
                continue

            linebuffer.append(line)
        
        logging.info(f'processed {current_chrom}')
        chrom_data = pd.DataFrame(
            linebuffer,
            columns = ['chrom', 'start', 'end', 'name']
        )
        yield chrom_data


def get_h5_group(h5, key_sequence):
    grp = h5
    for key in key_sequence:
        grp = grp[key]
    
    return grp


def write_annotation(grp, name, data, h5opts):
    if name in grp:
        del grp[name]
    
    grp.create_dataset(name, data=data, **h5opts)


def copy_and_annotate_cooler(source, dest, annotations, h5opts = None, mcoolfile = True):
    fileops.cp(source, dest)
    keys = ['/']
    ofile = dest

    if mcoolfile:
        ofile, keystring = dest.split('::')
        keys.extend(keystring[1:].split('/'))
        
    with h5py.File(ofile, 'r+') as h5:
        rootgrp = get_h5_group(h5, keys)
        grp = rootgrp['bins']
        h5opts = create._create._set_h5opts(h5opts)
        for col in annotations.columns:
            write_annotation(grp, col, annotations[col].values, h5opts)
    
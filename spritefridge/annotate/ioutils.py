from cooler import fileops, create

import h5py


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
    
import gc
import cooler
import h5py

import numpy as np

from scipy.sparse import csr_matrix, triu, tril
from cooler.util import parse_cooler_uri


def check_weight(cooleruri, weight_name):
    '''
    checks if weight_name already exist in cooler file

    :param cooleruri:   uri to a given cooleruri
    :param weight_name: name of the weight to check for

    :return:            True if weight already in cooler else False
    '''

    cool_path, group_path = parse_cooler_uri(cooleruri)
    weight_exists = False
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        if grp['bins'].get(weight_name):
            weight_exists = True

    return weight_exists


def store_weights(cooleruri, bias, weightname):
    '''
    stores an iterable of values as a new weight column in the given cooleruri
    with name set to wightname. code taken from cooler's cooler balance see also
    https://github.com/mirnylab/cooler/blob/master/cooler/cli/balance.py

    :param cooleruri:   uri to a given cooler
    :param bias:        iterable containing balancing weights for each genomic bin
    :param weightname:  name of the weight column

    :return:            None
    '''
    cool_path, group_path = parse_cooler_uri(cooleruri)
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        # add the bias column to the file
        h5opts = dict(compression='gzip', compression_opts=6)
        grp['bins'].create_dataset(weightname, data=bias, **h5opts)


def get_resolutons(coolerpath):
    '''
    returns all resolutions present in a MultiCooler file
    :param coolerpath:  path to MultiCooler file

    :return:            list of strings denoting the resolutions present
    '''
    with h5py.File(coolerpath, 'r+') as h5:
        return list(h5['resolutions'].keys())


def rename_weights(cooleruri, name_map):
    cool_path, group_path = parse_cooler_uri(cooleruri)
    with h5py.File(cool_path, 'r+') as h5:
        grp = h5[group_path]
        h5opts = dict(compression='gzip', compression_opts=6)
        for old_name, new_name in name_map.items():
            # add the bias column to the file
            weights = grp['bins'][old_name][()].copy()
            grp['bins'].create_dataset(new_name, data=weights, **h5opts)
            del grp['bins'][old_name]


def cooler_to_csr(cooleruri):
    '''
    loads a cooler into a csr matrix
    taken from HiCMatrix cool.py see also
    https://github.com/deeptools/HiCMatrix/blob/master/hicmatrix/lib/cool.py

    :param cooleruri:   uri to a given cooler

    :return:            data in cooler as scipy.sparse.csr_matrix
    '''
    cooler_file = cooler.Cooler(cooleruri)
    matrixDataFrame = cooler_file.matrix(
        balance=False, 
        sparse=True, 
        as_pixels=True
    )
    used_dtype = np.int32

    if np.iinfo(np.int32).max < cooler_file.info['nbins']:
        used_dtype = np.int64

    count_dtype = matrixDataFrame[0]['count'].dtype
    data = np.empty(
        cooler_file.info['nnz'], 
        dtype=count_dtype
    )
    instances = np.empty(
        cooler_file.info['nnz'], 
        dtype=used_dtype
    )
    features = np.empty(
        cooler_file.info['nnz'], 
        dtype=used_dtype
    )

    i = 0
    size = cooler_file.info['nbins'] // 32
    if size == 0:
        size = 1

    start_pos = 0
    while i < cooler_file.info['nbins']:
        matrixDataFrameChunk = matrixDataFrame[i:i + size]
        _data = matrixDataFrameChunk['count'].values.astype(count_dtype)
        _instances = matrixDataFrameChunk['bin1_id'].values.astype(used_dtype)
        _features = matrixDataFrameChunk['bin2_id'].values.astype(used_dtype)

        data[start_pos:start_pos + len(_data)] = _data
        instances[start_pos:start_pos + len(_instances)] = _instances
        features[start_pos:start_pos + len(_features)] = _features

        start_pos += len(_features)
        i += size

        del _data
        del _instances
        del _features

    matrix = csr_matrix(
        (data, (instances, features)),
        shape = (
            np.int(cooler_file.info['nbins']), 
            np.int(cooler_file.info['nbins'])
        ), 
        dtype = count_dtype
    )

    del data
    del instances
    del features
    gc.collect()

    # filling lower triangle in case only upper triangle was saved
    if tril(matrix, k=-1).sum() == 0:
        # this case means that the lower triangle of the
        # symmetric matrix (below the main diagonal)
        # is zero. In this case, replace the lower
        # triangle using the upper triangle
        matrix = matrix + triu(matrix, 1).T

    return matrix

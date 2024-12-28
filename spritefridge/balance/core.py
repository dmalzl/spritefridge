import logging
import cooler

import numpy as np
import multiprocessing as mp

from krbalancing import *
from .ioutils import cooler_to_csr
from cooler import ice


def remove_nan_bin_weights(csrmatrix, weights):
    '''
    removes weights from the balancing weight vector that belong to
    bins with 0 coverage (i.e. NaN bins)

    :param csrmatrix:   cooler as scipy.sparse.csr_matrix
    :param weights:     numpy.array of weights

    :return:            processed weight vector
    '''
    nnz_rows = csrmatrix.getnnz(1)
    nnz_cols = csrmatrix.getnnz(0)
    weights[(nnz_rows == 0) & (nnz_cols == 0)] = np.nan
    return weights


def balance_kr(cooleruri, per_chrom = False):
    '''
    applies KR matrix balancing to the given cooleruri
    code taken from HiCExplorer's hicCorrectMatrix see also
    https://github.com/deeptools/HiCExplorer/blob/master/hicexplorer/hicCorrectMatrix.py

    :param cooleruri:   uri to a given cooler

    :return:            KR balancing weights
    '''
    csrmatrix = cooler_to_csr(cooleruri)

    def apply_kr(m):
        kr = kr_balancing(
            m.shape[0],
            m.shape[1],
            m.count_nonzero(),
            m.indptr.astype(np.int64, copy=False),
            m.indices.astype(np.int64, copy=False),
            m.data.astype(np.float64, copy=False)
        )
        kr.computeKR()

        # set it to False since the vector is already normalised
        # with the previous True
        # correction_factors = np.true_divide(1, kr.get_normalisation_vector(False).todense())
        w = kr.get_normalisation_vector(False).todense()

        # flatten weights to comply with cooler format specification
        return np.array(w).flatten()

    if not per_chrom:
        weights = apply_kr(csrmatrix)

    else:
        coolermatrix = cooler.Cooler(cooleruri)
        weights = np.array([])
        for chrom in coolermatrix.chromnames:
            chromextent = coolermatrix.extent(chrom)
            chromcsr = csrmatrix[chromextent[0]: chromextent[1], chromextent[0]: chromextent[1]]
            chromweights = apply_kr(chromcsr)
            weights = np.concatenate([weights, np.array(chromweights).flatten()])

    return remove_nan_bin_weights(csrmatrix, weights)


def balance_ic(cooleruri, nproc, per_chrom):
    '''
    applies IC matrix balancing to a given cooleruri
    code taken from cooler's cooler balance see also
    https://github.com/mirnylab/cooler/blob/master/cooler/cli/balance.py
    :param cooleruri:   uri to a given cooler
    :param nproc:       number of processors to use for balancing

    :return:            IC balancing weights
    '''
    clr = cooler.Cooler(cooleruri)
    try:
        if nproc > 1:
            pool = mp.Pool(nproc)
            map_ = pool.imap_unordered
        else:
            map_ = map

        bias, stats = ice.iterative_correction(
            clr,
            chunksize=int(10e6),
            cis_only=per_chrom,
            trans_only=False,
            tol=1e-5,
            min_nnz=10,
            min_count=0,
            blacklist=None,
            mad_max=5,
            max_iters=500,
            ignore_diags=2,
            rescale_marginals=True,
            use_lock=False,
            map=map_
        )

    finally:
        if nproc > 1:
            pool.close()

    if not stats['converged']:
        logging.error('Iteration limit reached without convergence')
        logging.error('Storing final result. Check log to assess convergence.')

    return bias

import pandas as pd
import numpy as np
import itertools as it

from functools import partial


# all this could surely also be done by iterating over pixels
# however this needs additional overhead as pixels selector does not respect bin ids
# so one has to ensure that the pixel frames are aligned which is way easier with matrix


def merge_pixels(coolers, rowlo, rowhi, collo, colhi, round_floats):
    pixels = pd.DataFrame()
    tmp_count_cols = ['count_x', 'count_y']
    for size, cool in coolers.items():
        selector = cool.matrix(balance = False, as_pixels = True)
        tmp = selector[rowlo:rowhi, collo:colhi]
        if tmp.empty:
            continue

        tmp['count'] = tmp['count'].astype(float)
        tmp.loc[:, 'count'] = tmp['count'] * (2/size)

        if pixels.empty:
            pixels = tmp
            continue
        
        pixels = pixels.merge(
            tmp,
            on = ['bin1_id', 'bin2_id'],
            how = 'outer'
        )
        pixels['count'] = pixels[tmp_count_cols].sum(axis = 1)
        pixels.drop(
            columns = tmp_count_cols,
            inplace = True
        )
    
    # in case all pixel frames are empty
    if pixels.empty:
        return pixels

    if round_floats:
        pixels['float_count'] = pixels['count']
        pixels.loc[:, 'count'] = pixels.float_count.round(0)

    return pixels


def get_nbins(coolers):
    if isinstance(coolers, dict):
        key = np.random.choice(list(coolers.keys()))
        cool = coolers[key]

    cool = coolers
    return cool.info['nbins']


def chunked_pixels_iterator(pixels_selector, coolers, nchunks = 100, round_floats = True):
    chunkextents = np.linspace(
        0, 
        get_nbins(coolers), 
        nchunks, 
        dtype = int
    )
    for rowlo, rowhi in it.pairwise(chunkextents):
        for collo, colhi in it.pairwise(chunkextents):
            pixels = pixels_selector(
                coolers,
                rowlo, rowhi,
                collo, colhi,
                round_floats
            )
            if pixels.empty:
                continue

            yield pixels


def pixels_from_matrix_block(cooler, rowlo, rowhi, collo, colhi, round_floats):
    pixels = cooler.matrix(balance = False, as_pixels = True)[rowlo:rowhi, collo:colhi]

    if round_floats:
        pixels['float_count'] = pixels['count']
        pixels.loc[:, 'count'] = pixels.float_count.round(0)

    return pixels


chunked_pixels_merge = partial(
    chunked_pixels_iterator,
    merge_pixels
)

chunked_pixels_block = partial(
    chunked_pixels_iterator,
    pixels_from_matrix_block
)

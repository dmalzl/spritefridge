import glob
import re
import os
import resource
import logging

import pandas as pd

from cooler import Cooler


cs_regex = re.compile('_(?P<cs>[0-9]+)_')


def check_file_limit(n_files):
    softlimit, hardlimit = resource.getrlimit(
        resource.RLIMIT_NOFILE
    )

    required = n_files + 10
    if softlimit < n_files:
        logging.info(
            'open file limit is too small for current merger. ' +
            f'Is soft = {softlimit}, hard = {hardlimit}. ' +
            f'Needs to be at least {n_files}. Setting to {required} temporarily'
        )
        softlimit = required
        hardlimit = required

        resource.setrlimit(
            resource.RLIMIT_NOFILE,
            (softlimit, hardlimit)
        )


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

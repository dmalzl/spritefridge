import glob
import re
import os

import pandas as pd

from cooler import Cooler


cs_regex = re.compile('_(?P<cs>[0-9]+)_')


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

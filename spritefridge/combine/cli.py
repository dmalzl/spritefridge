def add_combine(subparser):
    parser = subparser.add_parser(
        'combine',
        help = '''
        merge multiple coolers generated from pairs of a given clustersize to a single one.
        This relies on the cooler names containing the clustersize.
        '''
    )
    parser.add_argument(
        '--input',
        '-i',
        help = 'path to directory containig clustersize coolers to merge. The name of the coolers must contain the clustersize _(?P<cs>[0-9]+)_',
        required = True
    )
    parser.add_argument(
        '--nchunks',
        help = 'the number of chunks to divide the data into when merging pixels. This is useful to finetune memory usage',
        default = 10,
        type = int
    )
    parser.add_argument(
        '--floatcounts',
        help = 'if set stores count column as float, else stores them in a separate column count then contains rounded float counts',
        action = 'store_true',
        default = False
    )
    parser.add_argument(
        '--outfile',
        '-o',
        help = 'path to file where the merged Cooler should be written to',
        required = True
    )

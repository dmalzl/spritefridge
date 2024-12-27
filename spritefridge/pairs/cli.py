def add_pairs(subparser):
    parser = subparser.add_parser(
        'pairs',
        help = '''
        generate pairs files for each cluster size from alignments. Alignments have to be filtered such that
        they only contain valid alignments and no multimappers
        '''
    )
    parser.add_argument(
        '--bam',
        '-b',
        help = 'BAM file containing aligned SPRITE-seq data',
        required = True
    )
    parser.add_argument(
        '--outprefix',
        '-o',
        help = 'prefix of the pairsfiles to write',
        required = True
    )
    parser.add_argument(
        '--clustersizelow',
        '-cl',
        help = 'minimum clustersize to consider',
        default = 2,
        type = int
    )
    parser.add_argument(
        '--clustersizehigh',
        '-ch',
        help = 'maximum clustersize to consider',
        default = 1000,
        type = int
    )
    parser.add_argument(
        '--separator',
        '-s',
        help = 'separator to use for extracting the barcode sequence from the readname',
        default = '['
    )
    
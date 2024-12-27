def add_extractbc(subparsers):
    parser = subparsers.add_parser(
        'extractbc',
        help = '''
        extract barcode sequences from raw SPRITE-seq reads and append them to the readname.
        by default this only writes reads having complete barcode sequences as defined by the given layouts
        '''
    )
    parser.add_argument(
        '--read1',
        '-r1',
        help = '(gzipped) fastq file containing sequence data for read1',
        required = True
    )
    parser.add_argument(
        '--read2',
        '-r2',
        help = '(gzipped) fastq file containing sequence data for read2',
        required = True
    )
    parser.add_argument(
        '--barcodes',
        '-bc',
        help = 'tab-separated file containing barcode information with columns category, bcname, bcseq',
        required = True
    )
    parser.add_argument(
        '--layout1',
        '-l1',
        help = 'barcode layout for read1 of the form category1|category2|...',
        required = True
    )
    parser.add_argument(
        '--layout2',
        '-l2',
        help = 'barcode layout for read2 of the form category1|category2|...',
        required = True
    )
    parser.add_argument(
        '--spacerlen',
        help = 'length of the spacer sequences if used',
        type = int,
        default = 6
    )
    parser.add_argument(
        '--laxity',
        help = 'number of bases to read into the current part of the read for matching barcodes',
        type = int,
        default = 6
    )
    parser.add_argument(
        '--mismatches',
        '-m',
        help = 'number of allowed mismatches per barcode category of the form category1:m1,category2:m2,...',
        required = True
    )
    parser.add_argument(
        '--output',
        '-o',
        help = 'file to write processed reads to',
        required = True
    )
    parser.add_argument(
        '--writefiltered',
        help = 'if set, writes reads with incomplete barcode set to a separate file',
        default = False,
        action = 'store_true'
    )
    parser.add_argument(
        '--writer2',
        help = 'if set, also writes r2 file which is usually only needed for barcode extraction',
        default = False,
        action = 'store_true'
    )
    parser.add_argument(
        '--processes',
        '-p',
        help = '''
        the number of processes to use for processing the reads. amounts to p - 2 extraction threads. e.g. 
        if p = 4 we have one main thread 2 extraction threads and 1 writer thread. if p = 1 no additional threads are spawned
        ''',
        default = 1,
        type = int
    )
    
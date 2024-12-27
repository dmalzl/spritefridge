import argparse as ap


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


def add_annotate(subparser):
    parser = subparser.add_parser(
        'annotate',
        help = '''
        annotate bins of one or multiple coolers with cluster information.
        i.e. if a read of a given cluster overlaps a given bin. annotated data is written to new file using old filename as prefix
        '''
    )    
    parser.add_argument(
        '--input',
        '-i',
        nargs = '+',
        help = 'one or more coolers (possibly mcool) to annotate with respective clusterinfo',
        required = True
    )
    parser.add_argument(
        '--bed',
        '-b',
        help = '''
        BEDfiles containing all valid SPRITE reads annotated with their cluster membership. 
        must be sorted by chrom, start
        ''',
        required = True
    )
    

def parse_args():
    parser = ap.ArgumentParser()
    sub = parser.add_subparsers(dest = 'subcommand')
    add_combine(sub)
    add_annotate(sub)
    add_extractbc(sub)
    add_pairs(sub)
    return parser.parse_args()


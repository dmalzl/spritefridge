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
        nargs = '+',
        help = '''
        one or more BEDfiles containing all valid SPRITE reads annotated with their cluster membership. 
        must be sorted by chrom, start. annotation columns will be named after bed basenames
        '''
    )
    
def add_balance(subparser):
    parser = subparser.add_parser(
        'balance',
        help = '''
        balance multi resolution coolers with KR and ICE algorithms
        '''
    )
    parser.add_argument(
        '-m', '--mcool',
        required = True,
        help = 'MultiCooler file to balance'
    )
    parser.add_argument(
        '-p', '--processors',
        default = 1,
        type = int,
        help = 'number of processors to use for IC balancing'
    )
    parser.add_argument(
        '--output',
        '-o',
        required = True,
        help = 'name of the output mcool file. input is copied there and appended with the weights'
    )
    parser.add_argument(
        '--overwrite',
        help = 'if set overwrites existing weight columns',
        default = False,
        action = 'store_true'
    )
    parser.add_argument(
        '--maxiter',
        default = 1000,
        help = 'maximum iterations for ICE algorithm'
    )

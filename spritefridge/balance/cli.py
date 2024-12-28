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

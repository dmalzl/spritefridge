from .cli import parse_args
from .extractbc import extract_barcodes
from .combine import combine_coolers
from .annotate import annotate_coolers
from .pairs import make_pairs

import logging


logging.basicConfig(
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)

cmds = {
    'combine': combine_coolers,
    'annotate': annotate_coolers,
    'extractbc': extract_barcodes,
    'pairs': make_pairs
}

def main():
    args = parse_args()
    cmds[args.subcommand](args)


if __name__ == '__main__':
    main()

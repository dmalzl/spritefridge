from .extractbc import extract_barcodes, add_extractbc
from .combine import combine_coolers, add_combine
from .annotate import annotate_coolers, add_annotate
from .pairs import make_pairs, add_pairs
from .balance import balance_mcool, add_balance

import logging

import argparse as ap


logging.basicConfig(
    format='%(asctime)s - %(message)s',
    level=logging.INFO
)


cmds = {
    'combine': combine_coolers,
    'annotate': annotate_coolers,
    'extractbc': extract_barcodes,
    'pairs': make_pairs,
    'balance': balance_mcool
}


def parse_args():
    parser = ap.ArgumentParser()
    sub = parser.add_subparsers(dest = 'subcommand')
    add_combine(sub)
    add_annotate(sub)
    add_extractbc(sub)
    add_pairs(sub)
    add_balance(sub)
    return parser.parse_args()


def main():
    args = parse_args()
    cmds[args.subcommand](args)


if __name__ == '__main__':
    main()

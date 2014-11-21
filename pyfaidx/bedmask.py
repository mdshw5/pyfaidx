import argparse
import sys
from pyfaidx import Fasta, bed_split


def mask_sequence(args):
    assert len(args.default_seq) == 1
    fasta = Fasta(args.fasta, mutable=True)
    for line in args.bed:
        rname, start, end = bed_split(line)
        if args.action == 'replace':
            fasta[rname][start:end] = (end - start) * args.default_seq
        elif args.action == 'lowercase':
            fasta[rname][start:end] = fasta[rname][start:end].lowercase()


def main():
    parser = argparse.ArgumentParser(description="Mask regions in BED file "
                                                 "from FASTA file")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('bed', type=argparse.FileType('r'), help="bed file of regions to mask")
    parser.add_argument('--default_seq', type=str, default='N',
                        help='default base for missing positions. default: %(default)s')
    parser.add_argument('-a', '--action', type=str, default='replace', choices=['replace', 'lowercase'],
                        help="masking action. default: %(default)s")
    args = parser.parse_args()
    mask_sequence(args)

if __name__ == "__main__":
    main()

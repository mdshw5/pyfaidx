#!/usr/bin/env python
import argparse
import sys
import os.path
import re
from pyfaidx import Faidx, Fasta, wrap_sequence, FetchError, ucsc_split, bed_split
from collections import Counter

keepcharacters = (' ', '.', '_')


def write_sequence(args):
    _, ext = os.path.splitext(args.fasta)
    if ext:
        ext = ext[1:]  # remove the dot from extension
    filt_function = re.compile(args.regex).search
    fasta = Fasta(args.fasta, default_seq=args.default_seq, strict_bounds=not args.lazy, split_char=args.delimiter, filt_function=filt_function)

    regions_to_fetch, split_function = split_regions(args)
    if not regions_to_fetch:
        regions_to_fetch = tuple(fasta.keys())

    header = False
    for region in regions_to_fetch:
        name, start, end = split_function(region)
        if args.split_files:  # open output file based on sequence name
            filename = '.'.join(str(e) for e in (name, start, end, ext) if e)
            filename = ''.join(c for c in filename if c.isalnum() or c in keepcharacters)
            outfile = open(filename, 'w')
        elif args.out:
            outfile = args.out
        else:
            outfile = sys.stdout
        try:
            if args.transform:
                if not header and args.transform == 'nucleotide':
                    outfile.write("name\tstart\tend\tA\tT\tC\tG\tN\n")
                    header = True
                outfile.write(transform_sequence(args, fasta, name, start, end))
            else:
                for line in fetch_sequence(args, fasta, name, start, end):
                    outfile.write(line)
        except FetchError as e:
            raise FetchError(e.msg.rstrip() + "Try setting --lazy.\n")
        if args.split_files:
            outfile.close()
    fasta.__exit__()


def fetch_sequence(args, fasta, name, start=None, end=None):
    line_len = fasta.faidx.index[name].lenc
    sequence = fasta[name][start:end]
    if args.complement:
        sequence = sequence.complement
    if args.reverse:
        sequence = sequence.reverse
    if args.no_names:
        pass
    elif args.full_names:
        yield ''.join(['>', fasta[name].long_name, '\n'])
    else:
        if start or end:
            yield ''.join(['>', sequence.longname, '\n'])
        else:
            yield ''.join(['>', sequence.name, '\n'])
    for line in wrap_sequence(line_len, sequence.seq):
        yield line


def mask_sequence(args):
    fasta = Fasta(args.fasta, mutable=True, split_char=args.delimiter)

    regions_to_fetch, split_function = split_regions(args)

    for region in regions_to_fetch:
        rname, start, end = split_function(region)
        if args.mask_with_default_seq:
            if start and end:
                span = end - start
            elif not start and not end:
                span = len(fasta[rname])
            else:
                span = len(fasta[rname][start:end])
            fasta[rname][start:end] = span * args.default_seq
        elif args.mask_by_case:
            fasta[rname][start:end] = fasta[rname][start:end].lowercase()


def split_regions(args):
    if args.bed:
        regions_to_fetch = args.bed
        split_function = bed_split
    else:
        regions_to_fetch = args.regions
        split_function = ucsc_split
    return (regions_to_fetch, split_function)


def transform_sequence(args, fasta, name, start=None, end=None):
    line_len = fasta.faidx.index[name].lenc
    s = fasta[name][start:end]
    if args.transform == 'bed':
        return '{name}\t{start}\t{end}\n'.format(name=s.name, start=s.start, end=s.end)
    elif args.transform == 'chromsizes':
        return '{name}\t{length}\n'.format(name=s.name, length=len(s))
    elif args.transform == 'nucleotide':
        nucs = Counter(dict([('A', 0), ('T', 0), ('C', 0), ('G', 0), ('N', 0)]))
        nucs.update(str(s).upper())
        return '{name}\t{start}\t{end}\t{A}\t{T}\t{C}\t{G}\t{N}\n'.format(name=s.name, start=s.start, end=s.end, **nucs)
    elif args.transform == 'transposed':
        return '{name}\t{start}\t{end}\t{seq}\n'.format(name=s.name, start=s.start, end=s.end, seq=str(s))



def main(ext_args=None):
    from pyfaidx import __version__
    parser = argparse.ArgumentParser(description="Fetch sequences from FASTA. If no regions are specified, all entries in the input file are returned. Input FASTA file must be consistently line-wrapped, and line wrapping of output is based on input line lengths.",
                                     epilog="Please cite: Shirley MD, Ma Z, Pedersen BS, Wheelan SJ. (2015) Efficient \"pythonic\" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196 https://dx.doi.org/10.7287/peerj.preprints.970v1")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('regions', type=str, nargs='*', help="space separated regions of sequence to fetch e.g. chr1:1-1000")
    parser.add_argument('-b', '--bed', type=argparse.FileType('r'), help="bed file of regions")
    parser.add_argument('-o', '--out', type=argparse.FileType('w'), help="output file name (default: stdout)")
    parser.add_argument('-i', '--transform', type=str, choices=('bed', 'chromsizes', 'nucleotide', 'transposed'), help="transform the requested regions into another format. default: %(default)s")
    parser.add_argument('-c', '--complement', action="store_true", default=False, help="complement the sequence. default: %(default)s")
    parser.add_argument('-r', '--reverse', action="store_true", default=False, help="reverse the sequence. default: %(default)s")
    names = parser.add_mutually_exclusive_group()
    names.add_argument('-n', '--no-names', action="store_true", default=False, help="omit sequence names from output. default: %(default)s")
    names.add_argument('-f', '--full-names', action="store_true", default=False, help="output full names including description. default: %(default)s")
    parser.add_argument('-x', '--split-files', action="store_true", default=False, help="write each region to a separate file (names are derived from regions)")
    parser.add_argument('-l', '--lazy', action="store_true", default=False, help="fill in --default-seq for missing ranges. default: %(default)s")
    parser.add_argument('-s', '--default-seq', type=check_seq_length, default='N', help='default base for missing positions and masking. default: %(default)s')
    parser.add_argument('-d', '--delimiter', type=str, default=None, help='delimiter for splitting names to multiple values (duplicate names will be discarded). default: %(default)s')
    parser.add_argument('-g', '--regex', type=str, default='.*', help='regular expression for filtering non-matching sequence names. default: %(default)s')
    masking = parser.add_mutually_exclusive_group()
    masking.add_argument('-m', '--mask-with-default-seq', action="store_true", default=False, help="mask the FASTA file using --default-seq default: %(default)s")
    masking.add_argument('-M', '--mask-by-case', action="store_true", default=False, help="mask the FASTA file by changing to lowercase. default: %(default)s")
    parser.add_argument('--version', action="version", version=__version__, help="print pyfaidx version number")
    # print help usage if no arguments are supplied
    if len(sys.argv)==1 and not ext_args:
        parser.print_help()
        sys.exit(1)
    elif ext_args:
        args = parser.parse_args(ext_args)
    else:
        args = parser.parse_args()

    if args.mask_with_default_seq or args.mask_by_case:
        mask_sequence(args)
    else:
        write_sequence(args)


def check_seq_length(value):
    if len(value) != 1:
        raise argparse.ArgumentTypeError("--default-seq value must be a single character!")
    return value

if __name__ == "__main__":
    main()

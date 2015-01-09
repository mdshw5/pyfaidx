#!/usr/bin/env python
"""Copyright (c) 2013 Matt Shirley

Permission is hereby granted, free of charge, to any person obtaining a copy of
this software and associated documentation files (the "Software"), to deal in
the Software without restriction, including without limitation the rights to
use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of
the Software, and to permit persons to whom the Software is furnished to do so,
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE."""

import argparse
import sys
import os.path
from pyfaidx import Faidx, Fasta, wrap_sequence, FetchError, ucsc_split, bed_split

keepcharacters = (' ', '.', '_')


def write_sequence(args):
    _, ext = os.path.splitext(args.fasta)
    if ext:
        ext = ext[1:]  # remove the dot from extension
    fasta = Fasta(args.fasta, default_seq=args.default_seq, strict_bounds=not args.lazy, split_char=args.delimiter)

    if args.bed:
        regions_to_fetch = args.bed
        split_function = bed_split
    else:
        regions_to_fetch = args.regions
        split_function = ucsc_split
        if not regions_to_fetch:
            regions_to_fetch = tuple(fasta.keys())

    for region in regions_to_fetch:
        name, start, end = split_function(region)
        if args.split_files:  # open output file based on sequence name
            filename = '.'.join(str(e) for e in (name, start + 1, end, ext) if e)
            filename = ''.join(c for c in filename if c.isalnum() or c in keepcharacters)
            outfile = open(filename, 'w')
        else:
            outfile = sys.stdout
        try:
            for line in fetch_sequence(args, fasta, name, start, end):
                outfile.write(line)
        except FetchError as e:
            sys.stderr.write("Error! {0} Try setting --lazy.\n".format(e.msg.rstrip()))
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
    if not args.no_names:
        if start or end:
            yield ''.join(['>', sequence.longname, '\n'])
        else:
            yield ''.join(['>', sequence.name, '\n'])
    for line in wrap_sequence(line_len, sequence.seq):
        yield line


def main():
    parser = argparse.ArgumentParser(description="Fetch sequences from FASTA. If no regions are specified, all entries in the input file are returned. Input FASTA file must be consistently line-wrapped, and line wrapping of output is based on input line lengths.")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('regions', type=str, nargs='*', help="space separated regions of sequence to fetch e.g. chr1:1-1000")
    parser.add_argument('-b', '--bed', type=argparse.FileType('r'), help="bed file of regions")
    parser.add_argument('--stats', action="store_true", default=False, help="print basic stats about the file. default: %(default)s")
    parser.add_argument('-c', '--complement', action="store_true", default=False, help="complement the sequence. default: %(default)s")
    parser.add_argument('-r', '--reverse', action="store_true", default=False, help="reverse the sequence. default: %(default)s")
    parser.add_argument('-n', '--no-names', action="store_true", default=False, help="print sequences without names. default: %(default)s")
    parser.add_argument('--split-files', action="store_true", default=False, help="write each region to a separate file (names are derived from regions)")
    parser.add_argument('--lazy', action="store_true", default=False, help="lazy region bounds checking - fill in default_seq for missing ranges. default: %(default)s")
    parser.add_argument('--default-seq', type=str, default='N', help='default base for missing positions. default: %(default)s')
    parser.add_argument('-d', '--delimiter', type=str, default=None, help='delimiter for splitting names to multiple values. default: %(default)s')
    # print help usage if no arguments are supplied
    if len(sys.argv)==1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    if args.stats:
        for key, value in Faidx(args.fasta).index.items():
            sys.stdout.write("{name}\t{length}\n".format(name=key, length=value['rlen']))
        sys.exit(0)

    write_sequence(args)

if __name__ == "__main__":
    main()

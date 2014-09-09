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
from pyfaidx import Fasta, wrap_sequence, FetchError


def write_sequence(args):
    fasta = Fasta(args.fasta, default_seq=args.default_seq, strict_bounds=not args.lazy)
    if args.bed:
        regions_to_fetch = args.bed
        split_function = bed_split
    else:
        regions_to_fetch = args.regions
        split_function = ucsc_split
    for region in regions_to_fetch:
        rname, start, end = split_function(region)
        try:
            for line in fetch_sequence(args, fasta, rname, start, end):
                sys.stdout.write(line)
        except FetchError as e:
            sys.stderr.write("Error! {0} Try setting --lazy.\n".format(e.msg.rstrip()))
    fasta.__exit__()


def bed_split(bed_entry):
    rname, start, end = bed_entry.rstrip().split()[:3]
    start, end = (int(start), int(end))
    return (rname, start, end)


def ucsc_split(region):
    region = region.split()[0]
    try:
        rname, interval = region.split(':')
    except ValueError:
        rname = region
        interval = None
    try:
        start, end = interval.split('-')
        start, end = (int(start) - 1, int(end))
    except (AttributeError, ValueError):
        start, end = (None, None)
    return (rname, start, end)


def fetch_sequence(args, fasta, rname, start=None, end=None):
    line_len = fasta.faidx.index[rname].lenc
    sequence = fasta[rname][start:end]
    if args.complement:
        sequence = sequence.complement
    if args.reverse:
        sequence = sequence.reverse
    if args.name:
        yield ''.join(['>', sequence.name, '\n'])
    for line in wrap_sequence(line_len, sequence.seq):
        yield line


def main():
    parser = argparse.ArgumentParser(description="Fetch sequence from "
                                                 "faidx-indexed FASTA")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('regions', type=str, nargs='?',
                        help="space separated regions of sequence to "
                        "fetch e.g. chr1:1-1000")
    parser.add_argument('-b', '--bed', type=argparse.FileType('r'), help="bed file of regions")
    parser.add_argument('-n', '--name', action="store_true", default=True,
                        help="print sequence names. default: %(default)s")
    parser.add_argument('--default_seq', type=str, default='N',
                        help='default base for missing positions. default: %(default)s')
    parser.add_argument('--lazy', action="store_true", default=False,
                        help="lazy region bounds checking - fill in default_seq for missing ranges. default: %(default)s")
    parser.add_argument('--complement', action="store_true", default=False,
                        help="comlement the sequence. default: %(default)s")
    parser.add_argument('--reverse', action="store_true", default=False,
                        help="reverse the sequence. default: %(default)s")
    args = parser.parse_args()
    write_sequence(args)

if __name__ == "__main__":
    main()

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
from pyfaidx import Fasta


def fetch(args):
    with Fasta(args.fasta) as fasta:
        regions = args.regions
        if args.list:
            with args.list as listfile:
                for region in listfile:
                    regions.append(region.rstrip())
        for region in regions:
            region = region.split()[0]
            try:
                rname, interval = region.split(':')
            except ValueError:
                rname = region
                interval = None
            try:
                start, end = interval.split('-')
                sequence = fasta[rname][int(start) - 1:int(end)]
            except (AttributeError, ValueError):
                sequence = fasta[rname][:]
            if args.name:
                sys.stdout.write(sequence.__repr__())
                sys.stdout.write('\n')
            else:
                sys.stdout.write(str(sequence))
                sys.stdout.write('\n')


def main():
    parser = argparse.ArgumentParser(description="Fetch sequence from "
                                                 "faidx-indexed FASTA")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('regions', type=str, nargs='*',
                        help="space separated regions of sequence to "
                        "fetch e.g. chr1:1-1000")
    parser.add_argument('-l', '--list', type=argparse.FileType('r'), help="list of regions, one per line")
    parser.add_argument('-n', '--name', action="store_true", default=True,
                        help="print sequence names")
    args = parser.parse_args()
    fetch(args)

if __name__ == "__main__":
    main()

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
from pyfaidx import *

def fetch(args):
    rname, interval = args.region.split(':')
    start, end = interval.split('-')
    with Genome(args.fasta) as genome:
        sequence = genome[rname][int(start) - 1:int(end)]
        if args.name:
            print(sequence)
        else:
            print(*sequence.seq, sep='')

def main():
    parser = argparse.ArgumentParser(description='Fetch sequence from faidx-indexed FASTA')
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('-r', '--region', type=str, help='region of sequence to fetch e.g. chr1:1-1000')
    parser.add_argument('-n', '--name', action="store_true", default=False, help='print sequence names')
    args = parser.parse_args()
    fetch(args)
    
if __name__ == "__main__":
    main()

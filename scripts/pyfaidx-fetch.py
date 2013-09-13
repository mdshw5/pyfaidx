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

import sys
import argparse
import re
import pyfaidx

def main(args):
    rname, interval = args.region.split(':')
    start, end = interval.split('-')
    kwargs = dict(zip(['rname', 'start', 'end'], [rname, int(start), int(end)]))
    with pyfaidx.faidx(args.fasta) as faidx:
        sequence = faidx.fetch(**kwargs)
        if args.name:
            sys.stdout.write('{name}\n{seq}\n'.format(name='>' + args.region, seq=sequence.seq))
        else:
            sys.stdout.write('{seq}\n'.format(seq=sequence.seq))

def parse_options():
    parser = argparse.ArgumentParser(description='Fetch sequence from faidx-indexed FASTA')
    parser.add_argument('fasta', type=str, help='faidx indexed FASTA file')
    parser.add_argument('-r', '--region', type=str, help='region of sequence to fetch e.g. chr1:1-1000')
    parser.add_argument('-n', '--name', action="store_true", default=False, help='print sequence names')
    args = parser.parse_args()
    return args
    
if __name__ == "__main__":
    args = parse_options()
    main(args)

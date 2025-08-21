#!/usr/bin/env python

import argparse
import sys
import os.path
import re
from pyfaidx import Fasta, wrap_sequence, FetchError, ucsc_split, bed_split, get_valid_filename
from collections import defaultdict
from operator import itemgetter
from heapq import nlargest
from itertools import repeat

def detect_fasta_newline(filepath):
    """Detect the newline style used in a FASTA file by reading the first non-header line."""
    with open(filepath, 'rb') as f:
        for line in f:
            if not line.startswith(b'>'):
                if line.endswith(b'\r\n'):
                    return '\r\n'
                elif line.endswith(b'\n'):
                    return '\n'
                elif line.endswith(b'\r'):
                    return '\r'
    return '\n'  # fallback

def write_sequence(args):
    _, ext = os.path.splitext(args.fasta)
    if ext:
        ext = ext[1:]  # remove the dot from extension

    filt_function = re.compile(args.regex).search

    if args.invert_match:
        filt_function = lambda x: not re.compile(args.regex).search(x)

    fasta = Fasta(args.fasta, default_seq=args.default_seq, key_function=eval(args.header_function), strict_bounds=not args.lazy, split_char=args.delimiter, filt_function=filt_function, read_long_names=args.long_names, rebuild=not args.no_rebuild)

    regions_to_fetch, split_function = split_regions(args)
    if not regions_to_fetch:
        regions_to_fetch = fasta.keys()

    header = False
    for region in regions_to_fetch:
        name, start, end = split_function(region)
        # allow the split_funtion to return None to signify input we should skip
        if name == None:
            continue
        if args.size_range:
            if start is not None and end is not None:
                sequence_len = end - start
            else:
                sequence_len = len(fasta[name])
            if args.size_range[0] > sequence_len or args.size_range[1] < sequence_len:
                continue
        if args.split_files:  # open output file based on sequence name
            filename = '.'.join(str(e) for e in (name, start, end, ext) if e)
            filename = get_valid_filename(filename)
            outfile = open(filename, 'w')
        elif args.out:
            outfile = args.out
        else:
            outfile = sys.stdout
        try:
            if args.transform:
                if not header and args.transform == 'nucleotide':
                    outfile.write("name\tstart\tend\tA\tT\tC\tG\tN\tothers\n")
                    header = True
                outfile.write(transform_sequence(args, fasta, name, start, end))
            else:
                for line in fetch_sequence(args, fasta, name, start, end):
                    outfile.write(line)
        except FetchError as e:
            raise FetchError(str(e) + " Try setting --lazy.\n")
        if args.split_files:
            outfile.close()
    fasta.__exit__()


def fetch_sequence(args, fasta, name, start=None, end=None):
    try:
        line_len = fasta.faidx.index[name].lenc
        if args.auto_strand and start > end and start is not None and end is not None:
            # flip (0, 1] coordinates
            sequence = fasta[name][end - 1:start + 1]
            sequence = sequence.reverse.complement
        else:
            sequence = fasta[name][start:end]
    except KeyError:
        sys.stderr.write("warning: {name} not found in file\n".format(**locals()))
        return
    if args.complement:
        sequence = sequence.complement
    if args.reverse:
        sequence = sequence.reverse
    if args.no_output:
        return
    if args.no_names:
        pass
    else:
        if (start or end) and not args.no_coords:
            yield ''.join(['>', sequence.fancy_name, '\n'])
        else:
            yield ''.join(['>', sequence.name, '\n'])
    newline = detect_fasta_newline(args.fasta)
    for line in wrap_sequence(line_len, sequence.seq, newline=newline):
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
            fasta[rname][start:end] = str(fasta[rname][start:end]).lower()


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
    if args.complement:
        s = s.complement
    if args.reverse:
        s = s.reverse
    if args.no_output:
        return
    if args.transform == 'bed':
        return '{name}\t{start}\t{end}\n'.format(name=s.name, start=s.start - 1 , end=s.end)
    elif args.transform == 'chromsizes':
        return '{name}\t{length}\n'.format(name=s.name, length=len(s))
    elif args.transform == 'nucleotide':
        ss = str(s).upper()
        nucs = defaultdict(int)
        nucs.update([(c, ss.count(c)) for c in set(ss)])
        A = nucs.pop('A', 0)
        T = nucs.pop('T', 0)
        C = nucs.pop('C', 0)
        G = nucs.pop('G', 0)
        N = nucs.pop('N', 0)
        # If there are other nucleotides, we will report them as well
        # Set others to 0 if no other nucleotides are present
        if not nucs:
            others = 0
        else:
            others = '|'.join([':'.join((k, str(v))) for k, v in nucs.items()])
        # Return the nucleotide counts in a tab-separated format
        return '{sname}\t{sstart}\t{send}\t{A}\t{T}\t{C}\t{G}\t{N}\t{others}\n'.format(sname=s.name, sstart=s.start, send=s.end, **locals())
    elif args.transform == 'transposed':
        return '{name}\t{start}\t{end}\t{seq}\n'.format(name=s.name, start=s.start, end=s.end, seq=str(s))



def main(ext_args=None):
    from pyfaidx import __version__
    parser = argparse.ArgumentParser(description="Fetch sequences from FASTA. If no regions are specified, all entries in the input file are returned. Input FASTA file must be consistently line-wrapped, and line wrapping of output is based on input line lengths.",
                                     epilog="Please cite: Shirley MD, Ma Z, Pedersen BS, Wheelan SJ. (2015) Efficient \"pythonic\" access to FASTA files using pyfaidx. PeerJ PrePrints 3:e1196 https://dx.doi.org/10.7287/peerj.preprints.970v1")
    parser.add_argument('fasta', type=str, help='FASTA file')
    parser.add_argument('regions', type=str, nargs='*', help="space separated regions of sequence to fetch e.g. chr1:1-1000")
    _input = parser.add_argument_group('input options')
    output = parser.add_argument_group('output options')
    header = parser.add_argument_group('header options')
    _input.add_argument('-b', '--bed', type=argparse.FileType('r'), help="bed file of regions (zero-based start coordinate)")
    output.add_argument('-o', '--out', type=argparse.FileType('w'), help="output file name (default: stdout)")
    output.add_argument('-i', '--transform', type=str, choices=('bed', 'chromsizes', 'nucleotide', 'transposed'), help="transform the requested regions into another format. default: %(default)s")
    output.add_argument('-c', '--complement', action="store_true", default=False, help="complement the sequence. default: %(default)s")
    output.add_argument('-r', '--reverse', action="store_true", default=False, help="reverse the sequence. default: %(default)s")
    output.add_argument('-y', '--auto-strand', action="store_true", default=False, help="reverse complement the sequence when start > end coordinate. default: %(default)s")
    output.add_argument('-a', '--size-range', type=parse_size_range, default=None, help='selected sequences are in the size range [low, high]. example: 1,1000 default: %(default)s')
    names = header.add_mutually_exclusive_group()
    names.add_argument('-n', '--no-names', action="store_true", default=False, help="omit sequence names from output. default: %(default)s")
    names.add_argument('-f', '--long-names', action="store_true", default=False, help="output full (long) names from the input fasta headers. default: headers are truncated after the first whitespace")
    header.add_argument('-t', '--no-coords', action="store_true", default=False, help="omit coordinates (e.g. chr:start-end) from output headers. default: %(default)s")
    output.add_argument('-x', '--split-files', action="store_true", default=False, help="write each region to a separate file (names are derived from regions)")
    output.add_argument('-l', '--lazy', action="store_true", default=False, help="fill in --default-seq for missing ranges. default: %(default)s")
    output.add_argument('-s', '--default-seq', type=check_seq_length, default=None, help='default base for missing positions and masking. default: %(default)s')
    header.add_argument('-d', '--delimiter', type=str, default=None, help='delimiter for splitting names to multiple values (duplicate names will be discarded). default: %(default)s')
    header.add_argument('-e', '--header-function', type=str, default='lambda x: x.split()[0]', help='python function to modify header lines e.g: "lambda x: x.split("|")[0]". default: %(default)s')
    header.add_argument('-u', '--duplicates-action', type=str, default="stop", choices=("stop", "first", "last", "longest", "shortest"), help='entry to take when duplicate sequence names are encountered. default: %(default)s')
    matcher = parser.add_argument_group('matching arguments')
    matcher.add_argument('-g', '--regex', type=str, default='.*', help='selected sequences are those matching regular expression. default: %(default)s')
    matcher.add_argument('-v', '--invert-match', action="store_true", default=False, help="selected sequences are those not matching 'regions' argument. default: %(default)s")
    masking = output.add_mutually_exclusive_group()
    masking.add_argument('-m', '--mask-with-default-seq', action="store_true", default=False, help="mask the FASTA file using --default-seq default: %(default)s")
    masking.add_argument('-M', '--mask-by-case', action="store_true", default=False, help="mask the FASTA file by changing to lowercase. default: %(default)s")
    output.add_argument('--no-output', action="store_true", default=False, help="do not output any sequence. default: %(default)s")
    parser.add_argument('--no-rebuild', action="store_true", default=False, help="do not rebuild the .fai index even if it is out of date. default: %(default)s")
    parser.add_argument('--version', action="version", version=__version__, help="print pyfaidx version number")
    # print help usage if no arguments are supplied
    if len(sys.argv)==1 and not ext_args:
        parser.print_help()
        sys.exit(1)
    elif ext_args:
        args = parser.parse_args(ext_args)
    else:
        args = parser.parse_args()
        
    if args.auto_strand:
        if args.complement:
            sys.stderr.write("--auto-strand and --complement are both set. Are you sure this is what you want?\n")
        if args.reverse:
            sys.stderr.write("--auto-strand and --reverse are both set. Are you sure this is what you want?\n")

    if args.mask_with_default_seq or args.mask_by_case:
        mask_sequence(args)
    else:
        write_sequence(args)

    

def check_seq_length(value):
    if value is None:
        pass # default value
    elif len(value) != 1:
        # user is passing a single character
        raise argparse.ArgumentTypeError("--default-seq value must be a single character!")
    return value

def parse_size_range(value):
    """ Size range argument should be in the form start,end and is end-inclusive. """
    if value is None:
        return value
    try:
        start, end = value.replace(' ', '').replace('\t', '').split(',')
    except (TypeError, ValueError, IndexError):
        raise ValueError
    return (int(start), int(end))


class Counter(dict):
    '''Dict subclass for counting hashable objects.  Sometimes called a bag
    or multiset.  Elements are stored as dictionary keys and their counts
    are stored as dictionary values.
    '''

    def __init__(self, iterable=None, **kwds):
        '''Create a new, empty Counter object.  And if given, count elements
        from an input iterable.  Or, initialize the count from another mapping
        of elements to their counts.
        '''
        self.update(iterable, **kwds)

    def __missing__(self, key):
        return 0


    def most_common(self, n=None):
        '''List the n most common elements and their counts from the most
        common to the least.  If n is None, then list all element counts.
        '''
        items = self.items()
        if n is None:
            return sorted(items, key=itemgetter(1), reverse=True)
        return nlargest(n, items, key=itemgetter(1))

    def elements(self):
        '''Iterator over elements repeating each as many times as its count.

        If an element's count has been set to zero or is a negative number,
        elements() will ignore it.

        '''
        for elem, count in self.items():
            for _ in repeat(None, count):
                yield elem

    # Override dict methods where the meaning changes for Counter objects.

    @classmethod
    def fromkeys(cls, iterable, v=None):
        raise NotImplementedError(
            'Counter.fromkeys() is undefined.  Use Counter(iterable) instead.')

    def update(self, iterable=None, **kwds):
        '''Like dict.update() but add counts instead of replacing them.

        Source can be an iterable, a dictionary, or another Counter instance.

        '''
        if iterable is not None:
            if hasattr(iterable, 'iteritems'):
                if self:
                    self_get = self.get
                    for elem, count in iterable.iteritems():
                        self[elem] = self_get(elem, 0) + count
                else:
                    dict.update(self, iterable) # fast path when counter is empty
            else:
                self_get = self.get
                for elem in iterable:
                    self[elem] = self_get(elem, 0) + 1
        if kwds:
            self.update(kwds)

    def copy(self):
        'Like dict.copy() but returns a Counter instance instead of a dict.'
        return Counter(self)

    def __delitem__(self, elem):
        'Like dict.__delitem__() but does not raise KeyError for missing values.'
        if elem in self:
            dict.__delitem__(self, elem)

    def __repr__(self):
        if not self:
            return '%s()' % self.__class__.__name__
        items = ', '.join(map('%r: %r'.__mod__, self.most_common()))
        return '%s({%s})' % (self.__class__.__name__, items)

    # Multiset-style mathematical operations discussed in:
    #       Knuth TAOCP Volume II section 4.6.3 exercise 19
    #       and at http://en.wikipedia.org/wiki/Multiset
    #
    # Outputs guaranteed to only include positive counts.
    #
    # To strip negative and zero counts, add-in an empty counter:
    #       c += Counter()

    def __add__(self, other):
        '''Add counts from two counters.

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] + other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __sub__(self, other):
        ''' Subtract count, but keep only results with positive counts.

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        result = Counter()
        for elem in set(self) | set(other):
            newcount = self[elem] - other[elem]
            if newcount > 0:
                result[elem] = newcount
        return result

    def __or__(self, other):
        '''Union is the maximum of value in either of the input counters.

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _max = max
        result = Counter()
        for elem in set(self) | set(other):
            newcount = _max(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

    def __and__(self, other):
        ''' Intersection is the minimum of corresponding counts.

        '''
        if not isinstance(other, Counter):
            return NotImplemented
        _min = min
        result = Counter()
        if len(self) < len(other):
            self, other = other, self
        for elem in filter(self.__contains__, other):
            newcount = _min(self[elem], other[elem])
            if newcount > 0:
                result[elem] = newcount
        return result

if __name__ == "__main__":
    main()

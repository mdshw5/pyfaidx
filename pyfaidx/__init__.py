"""
Fasta file -> Faidx -> Fasta -> FastaRecord -> Sequence
"""

from __future__ import division
import sys
import os
from six import PY2, PY3, string_types
from six.moves import zip_longest

if PY2:
    import string

class FastaIndexingError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class FetchError(Exception):
    def __init__(self, msg):
        self.msg = msg

    def __str__(self):
        return repr(self.msg)


class Sequence(object):
    """
    name = FASTA entry name
    seq = FASTA sequence
    start, end = coordinates of subsequence (optional)
    comp = boolean switch for complement property
    """
    def __init__(self, name='', seq='', start=None, end=None):
        self.name = name
        self.seq = seq
        self.start = start
        self.end = end
        self.comp = False
        assert isinstance(name, string_types)
        assert isinstance(seq, string_types)

    def __getitem__(self, n):
        if isinstance(n, slice):
            start, stop, step = n.indices(len(self))
            if stop == -1:
                stop = start
            else:
                stop = len(self) - stop
            if self.start and self.end:
                return self.__class__(self.name, self.seq[n.start:n.stop:n.step],
                                  self.start + start, self.end - stop)
            else:
                return self.__class__(self.name, self.seq[n.start:n.stop:n.step])
        elif isinstance(n, int):
            if n < 0:
                n = len(self) + n
            if self.start:
                return self.__class__(self.name, self.seq[n], self.start + n,
                                  self.start + n)
            else:
                return self.__class__(self.name, self.seq[n])

    def __str__(self):
        return self.seq

    def __neg__(self):
        """ Returns the reverse compliment of sequence
        >>> x = Sequence(name='chr1', seq='ATCGTA', start=1, end=6)
        >>> -x
        >chr1 (complement):6-1
        TACGAT
        """
        return self[::-1].complement

    def __repr__(self):
        if self.comp:
            name = '{rname} (complement)'.format(rname=self.name)
        else:
            name = self.name
        if self.start:
            return '\n'.join(['>{name}:{start}-{end}'.format(name=name,
                              start=self.start, end=self.end), self.seq])
        else:
            return '\n'.join(['>{name}'.format(name=name), self.seq])

    def __len__(self):
        """
        >>> len(Sequence('chr1', 'ACT'))
        3
        """
        return len(self.seq)

    @property
    def complement(self):
        """ Returns the compliment of self.
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> x.complement
        >chr1 (complement)
        TAGCAT
        """
        comp = self.__class__(self.name, complement(self.seq),
                              start=self.start, end=self.end)
        comp.comp = False if self.comp else True
        return comp

    @property
    def reverse(self):
        """ Returns the reverse of self.
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> x.reverse
        >chr1
        ATGCTA
        """
        return self[::-1]

    @property
    def gc(self):
        """ Return the GC content of seq as a float
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> y = round(x.gc, 2)
        >>> y == 0.33
        True
        """
        g = self.seq.count('G')
        g += self.seq.count('g')
        c = self.seq.count('C')
        c += self.seq.count('c')
        return (g + c) / len(self.seq)


class Faidx(object):
    """ A python implementation of samtools faidx FASTA indexing """
    def __init__(self, filename, key_function=None, as_raw=False):
        """
        filename: name of fasta file
        key_function: optional callback function which should return a unique key for the self.index dictionary when given rname.
        as_raw: optional parameter to specify whether to return sequences as a Sequence() object or as a raw string. Default: False (i.e. return a Sequence() object).
        """
        self.filename = filename
        self.file = open(filename, 'rb')
        self.indexname = filename + '.fai'
        self.key_function = key_function if key_function else lambda rname: rname
        self.as_raw = as_raw
        if os.path.exists(self.indexname):
            self.index = {}
            with open(self.indexname) as index:
                for line in index:
                    line = line.strip()
                    rname, rlen, offset, lenc, lenb = line.split('\t')
                    rname = self.key_function(rname)
                    if rname in self.index:
                        raise ValueError('Duplicate key "%s"'%rname)
                    self.index[rname] = {'rlen': int(rlen),
                                         'offset': int(offset),
                                         'lenc': int(lenc),
                                         'lenb': int(lenb)}

        else:
            self.build_fai(self.filename, self.indexname)
            self.__init__(filename, key_function)

    def __repr__(self):
        return 'Faidx("%s")' % (self.filename)

    @staticmethod
    def build_fai(filename, outfile):
        """ Build faidx index and write to ``outfile``.
        faidx index is in the format:
        rname\trlen\toffset\tlen\tblen """
        with open(outfile, 'w') as indexfile:
            with open(filename, 'r') as fastafile:
                rname = None  # reference sequence name
                offset = 0  # binary offset of end of current line
                rlen = 0  # reference character length
                blen = 0  # binary line length (includes newline)
                clen = 0  # character line length
                short_lines = []  # lines shorter than blen
                line_number = 0
                for line in fastafile:
                    line_blen = len(line)
                    line_clen = len(line.rstrip('\n\r'))
                    if (line[0] == '>') and (rname is None):
                        rname = line.rstrip()[1:].split()[0]
                        offset += line_blen
                        thisoffset = offset
                    elif line[0] != '>':
                        if blen == 0:
                            blen = line_blen
                        # only one short line should be allowed
                        # before we hit the next header
                        elif line_clen == 0:
                            sys.stderr.write("Warning: blank line in >{0} at "
                                             "line {1:n}.".format(rname,
                                                                  line_number +
                                                                  1))
                        elif blen > line_blen:
                            short_lines.append(line_number)
                            if len(short_lines) > 1:
                                indexfile.close()
                                os.remove(indexfile.name)
                                raise FastaIndexingError("Line length of fasta"
                                                         " file is not "
                                                         "consistent! "
                                    "Early short line found in >{0} at "
                                    "line {1:n}.".format(rname,
                                                         short_lines[0] + 1))
                        elif blen < line_blen:
                            raise FastaIndexingError("Line length of fasta "
                                                     "file is not consistent! "
                                                     "Long line found in >{0} "
                                                     "at line {1:n}."
                                                     "".format(rname,
                                                               line_number +
                                                               1))
                        offset += line_blen
                        if clen == 0:
                            clen = line_clen
                        rlen += line_clen
                    elif (line[0] == '>') and (rname is not None):
                        indexfile.write("{0}\t{1}\t"
                                        "{2}\t{3}\t{4}\n".format(rname,
                                                                 rlen,
                                                                 thisoffset,
                                                                 clen,
                                                                 blen))
                        blen = 0
                        rlen = 0
                        clen = 0
                        short_lines = []
                        rname = line.rstrip()[1:].split()[0]
                        offset += line_blen
                        thisoffset = offset
                    line_number += 1
                indexfile.write("{0}\t{1}\t{2}"
                                "\t{3}\t{4}\n".format(rname, rlen,
                                                      thisoffset, clen,
                                                      blen))

    def fetch(self, rname, start, end, strict_bounds=False):
        """ Fetch the sequence ``[start:end]`` from ``rname`` using 1-based coordinates
        1. Count newlines before start
        2. Count newlines to end
        3. Difference of 1 and 2 is number of newlines in [start:end]
        4. Seek to start position, taking newlines into account
        5. Read to end position, return sequence without newlines
        """
        try:
            entry = self.index[rname]
        except KeyError:
            raise FetchError("Requested rname {0} does not exist! "
                             "Please check your FASTA file.".format(rname))
        start = start - 1  # make coordinates [0,1)
        offset = entry.get('offset')
        rlen = entry.get('rlen')
        line_len = entry.get('lenc')
        line_blen = entry.get('lenb')
        seq_len = end - start
        newlines_total = int(rlen / line_len * (line_blen - line_len))
        newlines_before = int((start - 1) / line_len *
                              (line_blen - line_len)) if start > 0 else 0
        newlines_to_end = int(end / line_len * (line_blen - line_len))
        newlines_inside = newlines_to_end - newlines_before
        seq_blen = newlines_inside + seq_len
        bstart = offset + newlines_before + start
        bend = offset + newlines_total + rlen
        self.file.seek(bstart)
        if seq_blen < 0:
            return Sequence(name=rname, start=0, end=0)
        if bstart + seq_blen <= bend:
            s = self.file.read(seq_blen)
        elif strict_bounds:
            raise FetchError("Requested end coordinate is outside of {}. Set "
                             "strict_bounds=False to ignore.\n".format(rname))
        else:
            s = self.file.read(bend - bstart)
        seq = s.decode('utf-8')

        if self.as_raw:
            return seq.replace('\n', '')
        return Sequence(name=rname, start=int(start + 1),
                        end=int(end), seq=seq.replace('\n', ''))


    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()


class FastaRecord(object):
    def __init__(self, name, fa=None):
        self.name = name
        self._fa = fa

    def __getitem__(self, n):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        if isinstance(n, slice):
            start, stop, step = n.start, n.stop, n.step
            if not start:
                start = 0
            if not stop:
                stop = len(self)
            if stop < 0:
                stop = len(self) + stop
            if start < 0:
                start = len(self) + start
            return self._fa.get_seq(self.name, start + 1, stop)[::step]

        elif isinstance(n, int):
            if n < 0:
                n = len(self) + n
            return self._fa.get_seq(self.name, n + 1, n + 1)

    def __repr__(self):
        return 'FastaRecord("%s")' % (self.name)

    def __len__(self):
        """ Return length of chromosome """
        return self._fa.faidx.index[self.name]['rlen']

    def __str__(self):
        return str(self[:])


class Fasta(object):
    def __init__(self, filename, default_seq=None, key_function=None, as_raw=False, strict_bounds=False):
        """
        An object that provides a pygr compatible interface.
        filename: name of fasta file
        default_seq: if given, this base will always be returned if region is unavailable.
        key_function: optional callback function which should return a unique key for the self._records dictionary when given rname.
        as_raw: optional parameter to specify whether to return sequences as a Sequence() object or as a raw string. Default: False (i.e. return a Sequence() object).
        """
        self.filename = filename
        self.faidx = Faidx(filename, key_function=key_function, as_raw=as_raw)
        self._records = dict((rname, FastaRecord(rname, self)) for
                             rname in self.faidx.index.keys())
        self._default_seq = default_seq
        self.strict_bounds=strict_bounds

    def __contains__(self, record):
        """Return True if genome contains record."""
        return record in self._records

    def __getitem__(self, record):
        """Return a chromosome by its name."""
        if record not in self._records:
            self._records[record] = FastaRecord(record, self)
        return self._records[record]

    def __repr__(self):
        if self.faidx.as_raw:
            return 'Fasta("%s", as_raw=True)' % (self.filename)
        else:
            return 'Fasta("%s")' % (self.filename)

    def __iter__(self):
        for key in self.keys():
            yield self[key]

    def keys(self):
        return self._records.keys()

    def get_seq(self, chrom, start, end):
        """Return a sequence by record name and interval [start, end).

        Coordinates are 0-based, end-exclusive.
        """
        # Get sequence from real genome object and save result.
        return self.faidx.fetch(chrom, start, end, self.strict_bounds)

    def close(self):
        self.__exit__()

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.faidx.__exit__(*args)


def wrap_sequence(n, sequence, fillvalue=''):
    args = [iter(sequence)] * n
    for line in zip_longest(fillvalue=fillvalue, *args):
        yield ''.join(line + ("\n",))


def complement(seq):
    """ Returns the compliment of seq.
    >>> seq = 'ATCGTA'
    >>> complement(seq)
    'TAGCAT'
    """
    if PY3:
        table = str.maketrans('ACTGNactg', 'TGACNtgac')
    elif PY2:
        table = string.maketrans('ACTGNactg', 'TGACNtgac')
    return str(seq).translate(table)


def translate_chr_name(from_name, to_name):
    chr_name_map = dict(zip(from_name, to_name))

    def map_to_function(rname):
        return chr_name_map[rname]

    return map_to_function


if __name__ == "__main__":
    import doctest
    doctest.testmod()

"""
Fasta file -> Faidx -> Fasta -> FastaRecord -> Sequence
"""

from __future__ import division
import sys
import os
import itertools
from six import PY2, PY3, string_types

if PY2:
    import string
    
class FastaIndexingError(Exception):
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
            return self.__class__(self.name, self.seq[n.start:n.stop:n.step], self.start + start, self.end - stop)
        elif isinstance(n, int):
            if n < 0:
                n = len(self) + n
            return self.__class__(self.name, self.seq[n], self.start + n, self.start + n)
        
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
        if PY3:
            table = str.maketrans('ACTGN','TGACN')
        elif PY2:
            table = string.maketrans('ACTGN','TGACN')
        comp = self.__class__(self.name, str(self.seq).translate(table), start=self.start, end=self.end)
        comp.comp = False if self.comp else True
        return comp
    
    @property    
    def gc(self):
        """ Return the GC content of seq as a float
        >>> x = Sequence(name='chr1', seq='ATCGTA')
        >>> y = round(x.gc, 2)
        >>> y == 0.33
        True
        """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return (g + c) / len(self.seq)

class Faidx(object):
    """ A python implementation of samtools faidx FASTA indexing """
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, 'rb')
        self.indexname = filename + '.fai'
        if os.path.exists(self.indexname):
            self.index = {}
            with open(self.indexname) as index:
                for line in index:
                    line = line.strip()
                    rname, rlen, offset, lenc, lenb = line.split('\t')
                    self.index[rname] = {'rlen':int(rlen), 'offset':int(offset), 'lenc':int(lenc), 'lenb':int(lenb)}

        else:
            self.build_fai(self.filename, self.indexname)
            self.__init__(filename)
            
    def __repr__(self):
        return 'Faidx("%s")' % (self.filename)

    @staticmethod
    def build_fai(filename, outfile):
        """ Build faidx index and write to ``outfile``.
        faidx index is in the format:
        rname\trlen\toffset\tlen\tblen """
        with open(outfile, 'w') as indexfile:
            with open(filename, 'rb') as fastafile:
                rname = None ## reference sequence name
                offset = 0 ## binary offset of end of current line
                rlen = 0 ## reference character length
                blen = 0 ## binary line length (includes newline)
                clen = 0 ## character line length
                short_lines = [] ## lines shorter than blen
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
                        ## only one short line should be allowed
                        ## before we hit the next header
                        elif blen > line_blen:
                            short_lines.append(line_number)
                            if len(short_lines) > 1:
                                indexfile.close()
                                os.remove(indexfile.name)
                                raise FastaIndexingError("Line length of fasta file is not consistent! "
                                    "Early short line found in >{0} at line {1:n}.".format(rname, short_lines[0] + 1))
                        elif blen < line_blen:
                            raise FastaIndexingError("Line length of fasta file is not consistent! "
                                    "Long line found in >{0} at line {1:n}.".format(rname, line_number + 1))
                        offset += line_blen
                        if clen == 0:
                            clen = line_clen
                        rlen += line_clen
                    elif (line[0] == '>') and (rname is not None):
                        indexfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(rname, rlen, thisoffset, clen, blen))
                        blen = 0
                        rlen = 0
                        clen = 0
                        short_lines = []
                        rname = line.rstrip()[1:].split()[0]
                        offset += line_blen
                        thisoffset = offset
                    line_number += 1
                indexfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(rname, rlen, thisoffset, clen, blen))

    def fetch(self, rname, start, end):
        """ Fetch the sequence ``[start:end]`` from ``rname`` using 1-based coordinates
        1. Count newlines before start
        2. Count newlines to end
        3. Difference of 1 and 2 is number of newlines in [start:end]
        4. Seek to start position, taking newlines into account
        5. Read to end position, return sequence without newlines """
        try:
            entry = self.index[rname]
        except KeyError:
            sys.exit("Requested rname {0} does not exist! Please check your FASTA file.".format(rname))
        start = start - 1 ## make coordinates [0,1)
        offset = entry.get('offset')
        rlen = entry.get('rlen')
        line_len = entry.get('lenc')
        line_blen = entry.get('lenb')
        seq_len = end - start
        newlines_total = int(rlen / line_len * (line_blen - line_len))
        newlines_before = int((start - 1) / line_len * (line_blen - line_len)) if start > 0 else 0
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
        else:
            s = self.file.read(bend - bstart)
        seq = s.decode('utf-8')
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

class Fasta(object):
    def __init__(self, filename, default_seq=None):
        """
        An object that provides a pygr compatible interface.
        filename: fasta file
        default_seq: if given, this base will always be returned if
            region is unavailable.
        """
        self.filename = filename
        self.faidx = Faidx(filename)
        self._records = dict((rname, FastaRecord(rname, self)) for rname in self.faidx.index.keys())
        self._default_seq = default_seq

    def __contains__(self, record):
        """Return True if genome contains record."""
        return record in self._records

    def __getitem__(self, record):
        """Return a chromosome by its name."""
        if record not in self._records:
            self._records[record] = FastaRecord(record, self)
        return self._records[record]

    def __repr__(self):
        return 'Fasta("%s")' % (self.filename)
    
    def keys(self):
        return self._records.keys()

    def get_seq(self, chrom, start, end):
        """Return a sequence by record name and interval [start, end).

        Coordinates are 0-based, end-exclusive.
        """
        # Get sequence from real genome object and save result.
        return self.faidx.fetch(chrom, start, end)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.faidx.__exit__(*args)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

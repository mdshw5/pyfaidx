from __future__ import division
import sys
import mmap
import os
import itertools
from six import PY2, PY3, string_types

if PY2:
    import string

class Fasta(object):
    """ Hold name and sequence returned by `py:class:Reader` 
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
            return self.__class__(self.name, self.seq[n.start:n.stop])
        elif isinstance(n, key):
            return self.__class__(self.name, self.seq[n])
        
    def __str__(self):
        return self.seq
        
    def __neg__(self):
        """ Returns the reverse compliment of sequence 
        >>> x = Fasta(name='chr1', seq='ATCGTA', start=1, end=6)
        >>> -x
        >chr1 (complement)
        TACGAT
        """
        return self.__class__(self.name, self.seq[::-1], start=self.end, end=self.start).complement

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
        >>> len(Fasta('chr1', 'ACT'))
        3
        """
        return len(self.seq)
    
    @property
    def complement(self):
        """ Returns the compliment of self.
        >>> x = Fasta(name='chr1', seq='ATCGTA')
        >>> x.complement
        >chr1 (complement)
        TAGCAT
        """
        if PY3:
            table = str.maketrans('ACTGN','TGACN')
        elif PY2:
            table = string.maketrans('ACTGN','TGACN')
        comp = self.__class__(self.name, self.seq.translate(table), start=self.start, end=self.end)
        comp.comp = False if self.comp else True
        return comp

class Faidx(object):
    """ A python implementation of samtools faidx FASTA indexing """
    def __init__(self, filename):
        self.filename = filename
        self.file = open(filename, 'r+b')
        self.m = mmap.mmap(self.file.fileno(), 0)
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

    @staticmethod
    def build_fai(filename, outfile):
        """ Build faidx index and write to ``outfile``.
        faidx index is in the format:
        rname\trlen\toffset\tlen\tblen """
        with open(outfile, 'w') as indexfile:
            with open(filename, 'r') as fastafile:
                rname = None
                offset = 0
                rlen = 0
                blen = 0
                clen = 0
                for line in fastafile:
                    if (line[0] == '>') and (rname is None):
                        rname = line.rstrip()[1:].split()[0]
                        offset += len(line)
                        thisoffset = offset
                    elif line[0] != '>':
                        if blen == 0:
                            blen = len(line)
                        offset += len(line)
                        if clen == 0:
                            clen = len(line.rstrip())
                        rlen += len(line.rstrip())
                    elif (line[0] == '>') and (rname is not None):
                        indexfile.write("{0}\t{1}\t{2}\t{3}\t{4}\n".format(rname, rlen, thisoffset, clen, blen))
                        blen = 0
                        rlen = 0
                        rname = line.rstrip()[1:].split()[0]
                        offset += len(line)
                        thisoffset = offset
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
        self.m.seek(bstart)
        if seq_blen < 0:
            return Fasta(rname)
        if bstart + seq_blen <= bend:
            s = self.m.read(seq_blen)
        else:
            s = self.m.read(bend - bstart)
        seq = s.decode('utf-8')
        return Fasta(name=rname, start=int(start + 1),
            end=int(end), seq=seq.replace('\n', ''))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()

class Chromosome(object):
    def __init__(self, name, genome=None):
        self.name = name
        self.genome = genome

    def __getitem__(self, n):
        """Return sequence from region [start, end)

        Coordinates are 0-based, end-exclusive."""
        if isinstance(n, slice):
            return self.genome.get_seq(self.name, n.start + 1, n.stop)
        elif isinstance(n, key):
            return self.genome.get_seq(self.name, n.start + 1, n.start + 1)

    def __repr__(self):
        return 'Chromosome("%s")' % (self.name)


class Genome(object):
    def __init__(self, filename, default_seq=None):
        """
        A genome object that provides a pygr compatible interface.
        filename: fasta file
        default_seq: if given, this base will always be returned if
            region is unavailable.
        """
        self._genome = Faidx(filename)
        self._chroms = dict((rname, Chromosome(rname, self)) for rname in self._genome.index.keys())
        self._default_seq = default_seq

    def __contains__(self, chrom):
        """Return True if genome contains chromosome."""
        return chrom in self._chroms

    def __getitem__(self, chrom):
        """Return a chromosome by its name."""
        if chrom not in self._chroms:
            self._chroms[chrom] = Chromosome(chrom, self)
        return self._chroms[chrom]

    def get_seq(self, chrom, start, end):
        """Return a sequence by chromosome name and region [start, end).

        Coordinates are 0-based, end-exclusive.
        """
        # Get sequence from real genome object and save result.
        return self._genome.fetch(chrom, start, end)

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self._genome.__exit__()

def gc(seq):
    """ Return the GC content of seq as a float
    >>> x = 'ATCGTA'
    >>> y = round(gc(x), 2)
    >>> y == 0.33
    True
    """
    g = seq.count('G')
    c = seq.count('C')
    return (g + c) / len(seq)

if __name__ == "__main__":
    import doctest
    doctest.testmod()

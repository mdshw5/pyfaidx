from __future__ import division
import os
import sys
import gzip
import mmap

class reader:
    """ 
    A class to read the names and sequences from a fasta file.
    """
    def __init__(self, filename):
        name, ext = os.path.splitext(filename)
        if ext == '.gz':
            self.file = gzip.open(filename, 'rb')
        else:
            self.file = open(filename, 'r')

    def __iter__(self):
        """ 
        Return tuple (name, sequence).
        """
        seq = ''
        for line in self.file:
            if line[0] == '>':
                if seq != '':
                    yield fasta(name, seq.upper())
                    seq = ''
                name = line.rstrip()
            else:
                seq = seq + line.rstrip()
                
        yield fasta(name, seq.upper()) 

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()
        
class writer:
    """ Take a `py:class:fasta` object and file name, open file and write fasta """
    def __init__(self, filename=None, gz=False, stdout=False, wrap=70):
        if (gz == True) and (stdout == False):
            self.file = gzip.open(filename, 'wb')
        elif (gz == False) and (stdout == False):
            self.file = open(filename, 'w')
        elif stdout == True:
            self.file = sys.stdout
        elif (filename == None) and (stdout == False):
            raise IOError
        self.wrap = wrap

    def write(self, fasta):
        """ Write fastq with binary encoded conversion string appended to read strand """
        self.file.write(fasta.name + '\n')
        try:
            while len(fasta) > 0:
                self.file.write(fasta.seq[:self.wrap] + '\n')
                fasta.seq = fasta.seq[self.wrap:]
        except IndexError:
            self.file.write(fasta.seq + '\n')

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()
        
class fasta(object):
    """ Hold name and sequence returned by `py:class:fastaReader` """
    def __init__(self, name='', seq=''):
        self.name = name
        self.seq = seq
        
    def __getitem__(self, key):
        return self.__class__(self.name, self.seq[key])
        
    def __repr__(self):
        return '\n'.join([self.name, self.seq])
        
    def __len__(self):
        return len(self.seq)

    def reverse(self):
        """ Returns reverse ordered self """
        return self.__class__(self.name, self.seq[::-1])
        
    def complement(self):
        """ Returns the compliment of self. This only affects the sequence slot. """
        dict = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'N':'N'}
        compseq = ''.join(map(lambda x: dict[x], self.seq))
        return self.__class__(self.name, compseq)

    def revcomplement(self):
        """ Take the reverse compliment of self. """
        return self.reverse().complement()
        
    def gc(self):
        """ Return the GC content of self as a float """
        g = self.seq.count('G')
        c = self.seq.count('C')
        return (g + c) / len(self)
    
    def cpg(self):
        """ Return the number of CpG sites in self.seq """
        return self.seq.count('CG')

class faidx:
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
            self.build(self.filename, self.indexname)
            self.__init__(filename)

    def build(self, filename, outfile):
        """ Build faidx index and write to ``outfile``.
        faidx index is in the format:
        rname\trlen\toffset\tlen\tblen """
        with open(outfile, 'wb') as indexfile, open(filename, 'rb') as fastafile:
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
                    thisoffset = offset + len(line)
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
            print "rname does not exist"
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
            return fasta(rname, '')
        if bstart + seq_blen <= bend:
            s = self.m.read(seq_blen)
        else:
            s = self.m.read(bend - bstart)
        return fasta(rname, s.replace('\n',''))

    def __enter__(self):
        return self

    def __exit__(self, *args):
        self.file.close()
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()
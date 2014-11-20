import os
from pyfaidx import Fasta

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaRecordIter:
    def __init__(self):
        self.genes = os.path.join(path, 'data/genes.fasta')
        self.fasta = Fasta(self.genes)

    def test_fetch_whole_fasta(self):
        expect = open(self.genes).read()
        result = ''.join(['>' + record.name + '\n' + ''.join([line.seq + '\n' for line in record]) for record in self.fasta]).rstrip()
        assert expect == result

    def test_line_len(self):
        for record in self.fasta:
            assert len(next(iter(record))) == self.fasta.faidx.index[record.name].lenc

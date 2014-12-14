import os
from pyfaidx import Fasta

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaIntIndex:
    fasta = Fasta(os.path.join(path, 'data/genes.fasta'))

    def test_integer_index(self):
        expect = self.fasta['AB821309.1'][:100].seq
        result = self.fasta[0][:100].seq
        assert expect == result

import os
from pyfaidx import Faidx, Fasta
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureSequenceAsRaw:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.genes = Fasta(self.fasta)
        self.genes_as_raw = Fasta(self.fasta, as_raw=True)

    def test_as_raw_false(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = self.genes['KF435150.1'][100-1:150].seq.lower()
        assert result == expect

    def test_as_raw_true(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = self.genes_as_raw['KF435150.1'][100-1:150].lower()
        assert result == expect

    @raises(AttributeError)
    def test_as_raw_false_error(self):
        result = self.genes['KF435150.1'][100-1:150].lower()

    @raises(AttributeError)
    def test_as_raw_true_error(self):
        result = self.genes_as_raw['KF435150.1'][100-1:150].seq.lower()

    def test_as_raw_type_when_blen_lt_0(self):
        expect = ''
        result = self.genes_as_raw.faidx.fetch('KF435150.1', 10, 0)
        assert result == expect

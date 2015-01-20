import os
from pyfaidx import Faidx, Fasta, FetchError
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureBuffer:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.genes_buffer = Fasta(self.fasta, read_ahead=300, strict_bounds=True)
        self.genes = Fasta(self.fasta)

    def test_buffer_false(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = self.genes['KF435150.1'][100-1:150].seq.lower()
        assert result == expect

    def test_buffer_true(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = self.genes_buffer['KF435150.1'][100-1:150].seq.lower()
        assert result == expect

    def test_buffer_exceed(self):
        expect = 'atgacatcattttccacctctgctcagtgttcaacatctgacagtgcttgcaggatctctcctggacaaatcaatcaggtacgaccaaaactgccgcttttgaagattttgcatgcagcaggtgcgcaaggtgaaatgttcactgttaaagaggtcatgcactatttaggtcagtacataatggtgaagcaactttatgatcagcaggagcagcatatggtatattgtggtggagatcttttgggagaactactgggacgtcagagcttctccgtgaaagacccaagccctctctatgatatgctaagaaagaatcttgtcactttagccactgctactacagcaaagtgcagaggaaagttccacttccagaaaaagaactacagaagacgatatcccc'
        result = self.genes_buffer['KF435150.1'][0:400].seq.lower()
        assert result == expect

    @raises(FetchError)
    def test_bounds_error(self):
        result = self.genes_buffer['KF435150.1'][100-1:15000].seq.lower()

    @raises(ValueError)
    def test_buffer_value(self):
        Fasta(self.fasta, read_ahead = 0.5)

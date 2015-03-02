import os
from pyfaidx import Faidx, Fasta, FetchError
from nose.tools import raises
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureBuffer(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except IOError:
            pass  # some tests may delete this file

    def test_buffer_false(self):
        fasta = Fasta('data/genes.fasta', strict_bounds=True)
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].seq.lower()
        assert result == expect

    def test_buffer_true(self):
        fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
        result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].seq.lower()
        assert result == expect

    def test_buffer_exceed(self):
        fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
        expect = 'atgacatcattttccacctctgctcagtgttcaacatctgacagtgcttgcaggatctctcctggacaaatcaatcaggtacgaccaaaactgccgcttttgaagattttgcatgcagcaggtgcgcaaggtgaaatgttcactgttaaagaggtcatgcactatttaggtcagtacataatggtgaagcaactttatgatcagcaggagcagcatatggtatattgtggtggagatcttttgggagaactactgggacgtcagagcttctccgtgaaagacccaagccctctctatgatatgctaagaaagaatcttgtcactttagccactgctactacagcaaagtgcagaggaaagttccacttccagaaaaagaactacagaagacgatatcccc'
        result = fasta['gi|557361099|gb|KF435150.1|'][0:400].seq.lower()
        assert result == expect

    @raises(FetchError)
    def test_bounds_error(self):
        fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
        result = fasta['gi|557361099|gb|KF435150.1|'][100-1:15000].seq.lower()

    @raises(ValueError)
    def test_buffer_value(self):
        Fasta('data/genes.fasta', read_ahead=0.5)

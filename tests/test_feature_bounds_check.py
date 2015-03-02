import os
from pyfaidx import Faidx, FetchError
from nose.tools import raises
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureBoundsCheck:
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_fetch_whole_entry(self):
        faidx = Faidx('data/genes.fasta')
        expect = ('ATGACATCATTTTCCACCTCTGCTCAGTGTTCAACATCTGA'
                'CAGTGCTTGCAGGATCTCTCCTGGACAAATCAATCAGGTACGACCA'
                'AAACTGCCGCTTTTGAAGATTTTGCATGCAGCAGGTGCGCAAGG'
                'TGAAATGTTCACTGTTAAAGAGGTCATGCACTATTTAGGTCAGTACAT'
                'AATGGTGAAGCAACTTTATGATCAGCAGGAGCAGCATATGGTATATTG'
                'TGGTGGAGATCTTTTGGGAGAACTACTGGGACGTCAGAGCTTCTCCGTG'
                'AAAGACCCAAGCCCTCTCTATGATATGCTAAGAAAGAATCTTGTCACTTT'
                'AGCCACTGCTACTACAGCAAAGTGCAGAGGAAAGTTCCACTTCCAGAAAAA'
                'GAACTACAGAAGACGATATCCCCACACTGCCTACCTCAGAGCATAAATGCA'
                'TACATTCTAGAGAAGGTGATTGAAGTGGGAAAAAATGATGACCTGGAGGACTC')
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             1, 482)
        assert str(result) == expect

    def test_fetch_middle(self):
        faidx = Faidx('data/genes.fasta')
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             100, 150)
        assert str(result) == expect

    def test_fetch_end(self):
        faidx = Faidx('data/genes.fasta')
        expect = 'TC'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 482)
        assert str(result) == expect

    def test_fetch_border(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta')
        expect = 'TC'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 500)
        assert str(result) == expect

    def test_rev(self):
        faidx = Faidx('data/genes.fasta')
        expect = 'GA'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 482)
        assert str(-result) == expect, result

    @raises(FetchError)
    def test_fetch_past_bounds(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta', strict_bounds=True)
        expect = 'TC'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                                         480, 5000)

import os
from pyfaidx import Faidx, FetchError
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureBoundsCheck:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta)
        self.faidx_strict = Faidx(self.fasta, strict_bounds=True)

    def test_fetch_whole_entry(self):
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
        result = self.faidx.fetch('KF435150.1',
                             1, 482)
        assert str(result) == expect

    def test_fetch_middle(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.faidx.fetch('KF435150.1',
                             100, 150)
        assert str(result) == expect

    def test_fetch_end(self):
        expect = 'TC'
        result = self.faidx.fetch('KF435150.1',
                             480, 482)
        assert str(result) == expect

    def test_fetch_border(self):
        """ Fetch past the end of a gene entry """
        expect = 'TC'
        result = self.faidx.fetch('KF435150.1',
                             480, 500)
        assert str(result) == expect

    def test_rev(self):
        expect = 'GA'
        result = self.faidx.fetch('KF435150.1',
                             480, 482)
        assert str(-result) == expect, result

    @raises(FetchError)
    def test_fetch_past_bounds(self):
        """ Fetch past the end of a gene entry """
        expect = 'TC'
        result = self.faidx_strict.fetch('KF435150.1',
                                         480, 5000)

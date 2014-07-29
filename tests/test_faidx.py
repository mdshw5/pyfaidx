import os
from pyfaidx import *
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)

DEBUG = False

class TestFaidx:
    def __init__(self):
        self.fai = os.path.join(path, 'data/genes.fasta.fai')
        self.expect = os.path.join(path, 'data/expect/genes.fasta.fai')
        self.samtools = os.path.join(path, 'data/expect/genes.fasta.fai.samtools')
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta)

    def teardownclass(self):
        if DEBUG:
            return
        os.remove(self.fai)
        self.faidx.__exit__()

    def test_build(self):
        ### read results and assert they are the same as expected
        with open(self.fai, 'r') as fai:
            result = fai.read()
        with open(self.expect, 'r') as expect:
            expect = expect.read()
        assert result == expect

    def test_fetch_whole(self):
        """ Fetch the whole gene entry """
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
        """ Fetch the middle of a gene entry """
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.faidx.fetch('KF435150.1',
                             100, 150)
        assert str(result) == expect

    def test_fetch_end(self):
        """ Fetch the end of a gene entry """
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
        result = self.faidx.fetch('KF435150.1',
                             480, 5000, strict_bounds=True)

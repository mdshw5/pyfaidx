import os
from pyfaidx import Faidx, FastaIndexingError
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)


class TestIndexing:

    def __init__(self):
        self.fai = 'data/genes.fasta.fai'
        self.expect = 'data/expect/genes.fasta.fai'
        self.samtools = 'data/expect/genes.fasta.fai.samtools'
        self.fasta = 'data/genes.fasta'
        self.faidx = Faidx(self.fasta)

    def teardownclass(self):
        os.remove(self.faidx.indexname)

    def test_build(self):
        with open(self.expect, 'r') as fai:
            expect = fai.read()
        with open(self.faidx.indexname, 'r') as fai:
            result = fai.read()
        assert result == expect

    def test_samtools_compare(self):
        with open(self.samtools, 'r') as expect:
            expect = expect.read()
        with open(self.faidx.indexname, 'r') as fai:
            result = fai.read()
        assert result == expect

    def test_order(self):
        with open(self.expect, 'r') as fai:
            expect = [line.split()[0] for line in fai]
        result = list(self.faidx.index.keys())
        assert result == expect

    def test_remove_bad_index(self):
        try:
            Faidx('data/short_line.fasta')
        except FastaIndexingError:
            pass
        assert not os.path.exists('data/short_line.fasta.fai')

    @raises(FastaIndexingError)
    def test_short_line_lengths(self):
        Faidx('data/short_line.fasta')

    @raises(FastaIndexingError)
    def test_long_line_lengths(self):
        Faidx('data/long_line.fasta')

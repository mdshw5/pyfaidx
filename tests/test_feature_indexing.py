import os
from pyfaidx import Faidx

path = os.path.dirname(__file__)
os.chdir(path)

class TestIndexing:
    def __init__(self):
        self.fai = os.path.join(path, 'data/genes.fasta.fai')
        self.expect = os.path.join(path, 'data/expect/genes.fasta.fai')
        self.samtools = os.path.join(path, 'data/expect/genes.fasta.fai.samtools')
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta)

    def test_build(self):
        with open(self.expect, 'r') as fai:
            expect = fai.read()
        index = Faidx.build_fai(self.fasta)
        result = ''.join(index)
        assert result == expect

    def test_samtools_compare(self):
        with open(self.samtools, 'r') as expect:
            expect = expect.read()
        index = Faidx.build_fai(self.fasta)
        result = ''.join(index)
        assert result == expect

    def test_order(self):
        index = Faidx.build_fai(self.fasta)
        genes = [x.split()[0] for x in index]
        assert genes == list(self.faidx.index.keys())

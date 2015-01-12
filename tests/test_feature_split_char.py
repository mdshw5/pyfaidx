import os
from pyfaidx import Faidx, Fasta
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)


class TestFeatureSplitChar:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta, split_char='.')
        self.genes = Fasta(self.fasta, split_char='.')

    def test_keys(self):
        expect = ['3', 'AB821309', 'KF435149', 'KF435150', 'NM_000465', 'NM_001282543', 'NM_001282545', 'NM_001282548', 'NM_001282549', 'NR_104212', 'NR_104215', 'NR_104216', 'XM_005249642', 'XM_005249643', 'XM_005249644', 'XM_005249645', 'XM_005265507', 'XM_005265508', 'XR_241079', 'XR_241080', 'XR_241081']
        result = sorted(self.genes.keys())
        assert result == expect

    def test_key_function_by_dictionary_get_key(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.genes['KF435150'][100-1:150]
        assert str(result) == expect

    def test_key_function_by_fetch(self):
        expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
        result = self.faidx.fetch('KF435150',
                             100, 150)
        assert str(result) == expect

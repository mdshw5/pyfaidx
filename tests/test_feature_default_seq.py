import os
from pyfaidx import Faidx

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureDefaultSeq:
    def __init__(self):
        self.fasta = os.path.join(path, 'data/genes.fasta')
        self.faidx = Faidx(self.fasta, default_seq='N')

    def test_fetch_border_padded(self):
        """ Fetch past the end of a gene entry """
        expect = 'TCNNNNNNNNNNNNNNNNNNN'
        result = self.faidx.fetch('KF435150.1',
                             480, 500)
        assert str(result) == expect

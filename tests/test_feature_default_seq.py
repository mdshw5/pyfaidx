import os
from pyfaidx import Faidx
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFeatureDefaultSeq(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except FileNotFoundError:
            pass  # some tests may delete this file

    def test_fetch_border_padded(self):
        """ Fetch past the end of a gene entry """
        faidx = Faidx('data/genes.fasta', default_seq='N')
        expect = 'TCNNNNNNNNNNNNNNNNNNN'
        result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                             480, 500)
        assert str(result) == expect

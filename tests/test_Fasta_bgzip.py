import os
from pyfaidx import Fasta
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaBGZF(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.gz.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_integer_slice(self):
        fasta = Fasta('data/genes.fasta.gz')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][:100].seq
        result = fasta[0][:100].seq
        assert expect == result

    def test_integer_index(self):
        fasta = Fasta('data/genes.fasta.gz')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][100].seq
        result = fasta[0][100].seq
        assert expect == result

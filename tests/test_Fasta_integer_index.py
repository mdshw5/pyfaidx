import os
from pyfaidx import Fasta
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaIntIndex(TestCase):
    def setup_method(self):
        pass

    def teardown_method(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_integer_slice(self):
        fasta = Fasta('data/genes.fasta')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][:100].seq
        result = fasta[0][:100].seq
        assert expect == result

    def test_integer_index(self):
        fasta = Fasta('data/genes.fasta')
        expect = fasta['gi|563317589|dbj|AB821309.1|'][100].seq
        result = fasta[0][100].seq
        assert expect == result

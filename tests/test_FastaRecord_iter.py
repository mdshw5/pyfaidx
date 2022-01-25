import os
from pyfaidx import Fasta
from itertools import chain
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaRecordIter(TestCase):
    def setup_method(self):
        pass

    def teardown_method(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_fetch_whole_fasta(self):
        expect = [line.rstrip('\n') for line in open('data/genes.fasta') if line[0] != '>']
        result = list(chain(*([line for line in record] for record in Fasta('data/genes.fasta', as_raw=True))))
        assert expect == result

    def test_line_len(self):
        fasta = Fasta('data/genes.fasta')
        for record in fasta:
            assert len(next(iter(record))) == fasta.faidx.index[record.name].lenc

    def test_reverse_iter(self):
        expect = list(chain(*([line[::-1] for line in record][::-1] for record in Fasta('data/genes.fasta', as_raw=True))))
        result = list(chain(*([line for line in reversed(record)] for record in Fasta('data/genes.fasta', as_raw=True))))
        for a, b in zip(expect, result):
            print(a, b)
        assert expect == result

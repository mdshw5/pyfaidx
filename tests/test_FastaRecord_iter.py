import os
from pyfaidx import Fasta
from itertools import chain
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestFastaRecordIter(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except IOError:
            pass  # some tests may delete this file

    def test_fetch_whole_fasta(self):
        expect = [line.rstrip('\n') for line in open('data/genes.fasta') if line[0] != '>']
        result = list(chain(*([line for line in record] for record in Fasta('data/genes.fasta', as_raw=True))))
        assert expect == result

    def test_line_len(self):
        fasta = Fasta('data/genes.fasta')
        for record in fasta:
            assert len(next(iter(record))) == fasta.faidx.index[record.name].lenc

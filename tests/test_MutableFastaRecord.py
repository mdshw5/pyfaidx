import os
from pyfaidx import Fasta
from tempfile import NamedTemporaryFile
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestMutableFastaRecord(TestCase):
    def setUp(self):
        with open('data/genes_mutable.fasta', 'wb') as mutable:
            mutable.write(open('data/genes.fasta', 'rb').read())
        self.mutable_fasta = Fasta('data/genes_mutable.fasta', mutable=True)

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except IOError:
            pass  # some tests may delete this file
        try:
            os.remove('data/genes_mutable.fasta')
        except IOError:
            pass  # some tests may delete this file
        try:
            os.remove('data/genes_mutable.fasta.fai')
        except IOError:
            pass  # some tests may delete this file

    def test_mutate_fasta_to_same(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        fasta = Fasta('data/genes.fasta', mutable=False)
        chunk = fasta['gi|557361099|gb|KF435150.1|'][0:100]
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk.seq
        assert str(fasta['gi|557361099|gb|KF435150.1|']) == str(mutable['gi|557361099|gb|KF435150.1|'])

    def test_mutate_fasta_to_N(self):
        mutable = Fasta('data/genes_mutable.fasta', mutable=True)
        chunk = 100 * 'N'
        mutable['gi|557361099|gb|KF435150.1|'][0:100] = chunk
        assert mutable['gi|557361099|gb|KF435150.1|'][0:100].seq == chunk

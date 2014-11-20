import os
from pyfaidx import Fasta
from tempfile import NamedTemporaryFile

path = os.path.dirname(__file__)
os.chdir(path)

class TestMutableFastaRecord:
    def __init__(self):
        self.genes = os.path.join(path, 'data/genes.fasta')
        self.fasta = Fasta(self.genes)

    def setup(self):
        self.genes_copy = NamedTemporaryFile(mode='rb+')
        self.genes_copy.write(open(self.genes, 'rb').read())
        self.genes_copy.seek(0)
        self.mutable_fasta = Fasta(self.genes_copy.name, mutable=True)

    def teardown(self):
        self.genes_copy.close()  # deletes temporary file

    def test_mutate_fasta(self):
        chunk = self.fasta['KF435150.1'][0:100]
        self.mutable_fasta['KF435150.1'][0:100] = chunk.seq
        assert str(self.fasta['KF435150.1']) == str(self.mutable_fasta['KF435150.1'])

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
        self.genes_copy = NamedTemporaryFile(mode='wb', delete=False)
        self.genes_copy.write(open(self.genes, 'rb').read())
        self.genes_copy.close()
        self.mutable_fasta = Fasta(self.genes_copy.name, mutable=True)

    def teardown(self):
        self.mutable_fasta.__exit__()
        os.remove(self.mutable_fasta.filename)  # deletes temporary file

    def test_mutate_fasta_to_same(self):
        chunk = self.fasta['KF435150.1'][0:100]
        self.mutable_fasta['KF435150.1'][0:100] = chunk.seq
        assert str(self.fasta['KF435150.1']) == str(self.mutable_fasta['KF435150.1'])

    def test_mutate_fasta_to_N(self):
        chunk = 100 * 'N'
        self.mutable_fasta['KF435150.1'][0:100] = chunk
        assert self.mutable_fasta['KF435150.1'][0:100].seq == chunk

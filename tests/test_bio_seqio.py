import os
from pyfaidx import Fasta, FetchError
from Bio import SeqIO

path = os.path.dirname(__file__)
os.chdir(path)

class TestBioSeqIO:
    def __init__(self):
        self.genes = os.path.join(path, 'data/genes.fasta')
        self.fasta = Fasta(self.genes)
        with open(self.genes, "rU") as fh:
            self.seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

    def test_fetch_whole_entry(self):
        assert str(self.fasta['KF435150.1']) == str(self.seqio['KF435150.1'].seq)
        assert self.fasta['KF435150.1'].name == str(self.seqio['KF435150.1'].name)

    def test_slice_whole_entry(self):
        assert str(self.fasta['KF435150.1'][::3]) == str(self.seqio['KF435150.1'].seq[::3])

    def test_revcomp_whole_entry(self):
        assert str(self.fasta['KF435150.1'][:].reverse.complement) == str(self.seqio['KF435150.1'].reverse_complement().seq)

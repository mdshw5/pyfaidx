import os
from pyfaidx import Fasta, FetchError
from nose.plugins.skip import Skip, SkipTest
try:
    from Bio import SeqIO
    test_bio = True
except ImportError:
    test_bio = False

path = os.path.dirname(__file__)
os.chdir(path)

class TestBioSeqIO:
    def __init__(self):
        self.genes = os.path.join(path, 'data/genes.fasta')
        self.fasta = Fasta(self.genes)
        if test_bio:
            with open(self.genes, "rU") as fh:
                self.seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))

    def test_fetch_whole_entry(self):
        if test_bio:
            assert str(self.fasta['KF435150.1']) == str(self.seqio['KF435150.1'].seq)
            assert self.fasta['KF435150.1'].name == str(self.seqio['KF435150.1'].name)
        else:
            raise SkipTest

    def test_slice_whole_entry(self):
        if test_bio:
            assert str(self.fasta['KF435150.1'][::3]) == str(self.seqio['KF435150.1'].seq[::3])
        else:
            raise SkipTest

    def test_revcomp_whole_entry(self):
        if test_bio:
            assert str(self.fasta['KF435150.1'][:].reverse.complement) == str(self.seqio['KF435150.1'].reverse_complement().seq)
        else:
            raise SkipTest

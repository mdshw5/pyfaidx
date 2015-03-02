import os
from pyfaidx import Fasta, FetchError
from nose.plugins.skip import Skip, SkipTest
from unittest import TestCase
try:
    from Bio import SeqIO
    test_bio = True
except ImportError:
    test_bio = False

path = os.path.dirname(__file__)
os.chdir(path)

class TestBioSeqIO(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except IOError:
            pass  # some tests may delete this file

    def test_fetch_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        if test_bio:
            with open('data/genes.fasta', "rU") as fh:
                seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            assert str(fasta['gi|557361099|gb|KF435150.1|']) == str(seqio['gi|557361099|gb|KF435150.1|'].seq)
            assert fasta['gi|557361099|gb|KF435150.1|'].name == str(seqio['gi|557361099|gb|KF435150.1|'].name)
        else:
            raise SkipTest

    def test_slice_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        if test_bio:
            with open('data/genes.fasta', "rU") as fh:
                seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            assert str(fasta['gi|557361099|gb|KF435150.1|'][::3]) == str(seqio['gi|557361099|gb|KF435150.1|'].seq[::3])
        else:
            raise SkipTest

    def test_revcomp_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        if test_bio:
            with open('data/genes.fasta', "rU") as fh:
                seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
            assert str(fasta['gi|557361099|gb|KF435150.1|'][:].reverse.complement) == str(seqio['gi|557361099|gb|KF435150.1|'].reverse_complement().seq)
        else:
            raise SkipTest

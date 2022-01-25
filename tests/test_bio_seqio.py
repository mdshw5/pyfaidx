import os
import pytest
from pyfaidx import Fasta, FetchError
from unittest import TestCase

path = os.path.dirname(__file__)
os.chdir(path)

class TestBioSeqIO(TestCase):
    def setUp(self):
        try:
            from Bio import SeqIO
        except ImportError:
            pytest.skip("biopython not installed.")

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_fetch_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        with open('data/genes.fasta', "rU") as fh:
            seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        assert str(fasta['gi|557361099|gb|KF435150.1|']) == str(seqio['gi|557361099|gb|KF435150.1|'].seq)
        assert fasta['gi|557361099|gb|KF435150.1|'].name == str(seqio['gi|557361099|gb|KF435150.1|'].name)

    def test_slice_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        with open('data/genes.fasta', "rU") as fh:
            seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        assert str(fasta['gi|557361099|gb|KF435150.1|'][::3]) == str(seqio['gi|557361099|gb|KF435150.1|'].seq[::3])

    def test_revcomp_whole_entry(self):
        fasta = Fasta('data/genes.fasta')
        with open('data/genes.fasta', "rU") as fh:
            seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
        assert str(fasta['gi|557361099|gb|KF435150.1|'][:].reverse.complement) == str(seqio['gi|557361099|gb|KF435150.1|'].reverse_complement().seq)
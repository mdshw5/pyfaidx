import os
import pytest
from pyfaidx import Fasta, FetchError

path = os.path.dirname(__file__)
os.chdir(path)

try:
    from Bio import SeqIO
    bio = True
except ImportError:
    bio = False
    
@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

@pytest.mark.skipif(not bio, reason="Biopython is not installed.")
def test_fetch_whole_entry(remove_index):
    fasta = Fasta('data/genes.fasta')
    with open('data/genes.fasta', "r") as fh:
        seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
    assert str(fasta['KF435150.1']) == str(seqio['KF435150.1'].seq)
    assert fasta['KF435150.1'].name == str(seqio['KF435150.1'].name)

@pytest.mark.skipif(not bio, reason="Biopython is not installed.")
def test_slice_whole_entry(remove_index):
    fasta = Fasta('data/genes.fasta')
    with open('data/genes.fasta', "r") as fh:
        seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
    assert str(fasta['KF435150.1'][::3]) == str(seqio['KF435150.1'].seq[::3])

@pytest.mark.skipif(not bio, reason="Biopython is not installed.")
def test_revcomp_whole_entry(remove_index):
    fasta = Fasta('data/genes.fasta')
    with open('data/genes.fasta', "r") as fh:
        seqio = SeqIO.to_dict(SeqIO.parse(fh, "fasta"))
    assert str(fasta['KF435150.1'][:].reverse.complement) == str(seqio['KF435150.1'].reverse_complement().seq)
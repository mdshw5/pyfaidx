import os
import pytest
from pyfaidx import Faidx, Fasta

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_as_raw_false(remove_index):
    fasta = Fasta('data/genes.fasta')
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
    result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].seq.lower()
    assert result == expect

def test_as_raw_true(remove_index):
    fasta = Fasta('data/genes.fasta', as_raw=True)
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
    result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].lower()
    assert result == expect

def test_as_raw_false_error(remove_index):
    fasta = Fasta('data/genes.fasta')
    with pytest.raises(AttributeError):
        result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].lower()

def test_as_raw_true_error(remove_index):
    fasta = Fasta('data/genes.fasta', as_raw=True)
    with pytest.raises(AttributeError):
        result = fasta['gi|557361099|gb|KF435150.1|'][100-1:150].seq.lower()

def test_as_raw_type_when_blen_lt_0(remove_index):
    fasta = Fasta('data/genes.fasta', as_raw=True)
    expect = ''
    result = fasta.faidx.fetch('gi|557361099|gb|KF435150.1|', 10, 0)
    assert result == expect

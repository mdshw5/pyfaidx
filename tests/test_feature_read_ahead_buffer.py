import os
import pytest
from pyfaidx import Faidx, Fasta, FetchError

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
        os.remove('data/genes.fasta.lower.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_buffer_false(remove_index):
    fasta = Fasta('data/genes.fasta', strict_bounds=True)
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
    result = fasta['KF435150.1'][100-1:150].seq.lower()
    assert result == expect

def test_buffer_true(remove_index):
    fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'.lower()
    result = fasta['KF435150.1'][100-1:150].seq.lower()
    assert result == expect

def test_buffer_exceed(remove_index):
    fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
    expect = 'atgacatcattttccacctctgctcagtgttcaacatctgacagtgcttgcaggatctctcctggacaaatcaatcaggtacgaccaaaactgccgcttttgaagattttgcatgcagcaggtgcgcaaggtgaaatgttcactgttaaagaggtcatgcactatttaggtcagtacataatggtgaagcaactttatgatcagcaggagcagcatatggtatattgtggtggagatcttttgggagaactactgggacgtcagagcttctccgtgaaagacccaagccctctctatgatatgctaagaaagaatcttgtcactttagccactgctactacagcaaagtgcagaggaaagttccacttccagaaaaagaactacagaagacgatatcccc'
    result = fasta['KF435150.1'][0:400].seq.lower()
    assert result == expect

def test_bounds_error(remove_index):
    fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
    with pytest.raises(FetchError):
        result = fasta['KF435150.1'][100-1:15000].seq.lower()

def test_buffer_value(remove_index):
    with pytest.raises(ValueError):
        Fasta('data/genes.fasta', read_ahead=0.5)

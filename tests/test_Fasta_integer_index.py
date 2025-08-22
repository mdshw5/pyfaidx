import os
import pytest
from pyfaidx import Fasta

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file
    
def test_integer_slice():
    fasta = Fasta('data/genes.fasta')
    expect = fasta['AB821309.1'][:100].seq
    result = fasta[0][:100].seq
    assert expect == result

def test_integer_index(remove_index):
    fasta = Fasta('data/genes.fasta')
    expect = fasta['AB821309.1'][100].seq
    result = fasta[0][100].seq
    assert expect == result

@pytest.mark.xfail(raises=IndexError)
def test_integer_index_out_of_bounds(remove_index):
    fasta = Fasta('data/genes.fasta')
    # This should raise an error because the index is out of bounds
    result = fasta[100]

@pytest.mark.xfail(raises=TypeError)
def test_integer_index_invalid_type(remove_index):
    fasta = Fasta('data/genes.fasta')
    # This should raise an error because the index is not an integer or string
    result = fasta[0.5]

@pytest.mark.xfail(raises=KeyError)
def test_integer_index_invalid_key(remove_index):
    fasta = Fasta('data/genes.fasta')
    # This should raise an error because the key does not exist
    result = fasta['non_existent_key']
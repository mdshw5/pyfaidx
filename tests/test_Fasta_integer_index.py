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
    
def test_integer_slice(remove_index):
    fasta = Fasta('data/genes.fasta')
    expect = fasta['gi|563317589|dbj|AB821309.1|'][:100].seq
    result = fasta[0][:100].seq
    assert expect == result

def test_integer_index(remove_index):
    fasta = Fasta('data/genes.fasta')
    expect = fasta['gi|563317589|dbj|AB821309.1|'][100].seq
    result = fasta[0][100].seq
    assert expect == result

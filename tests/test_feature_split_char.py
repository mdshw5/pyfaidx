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
    
def test_keys(remove_index):
    fasta = Fasta('data/genes.fasta', split_char='.', duplicate_action="drop")
    expect = ['AB821309', 'KF435149', 'KF435150', 'NR_104158', 'NR_104160', 'NR_104205', 'NR_104206', 'NR_104208', 'NR_104214', 'NR_104215', 'NR_104216', 'XR_241053', 'XR_241054', 'XR_241055', 'XR_241064', 'XR_241065', 'XR_241078', 'XR_241079', 'XR_241080', 'XR_241081']
    result = sorted(fasta.keys())
    assert result == expect

def test_key_function_by_dictionary_get_key(remove_index):
    fasta = Fasta('data/genes.fasta', split_char='.', duplicate_action="drop")
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = fasta['KF435150'][100-1:150]
    assert str(result) == expect

def test_key_function_by_fetch(remove_index):
    faidx = Faidx('data/genes.fasta', split_char='.', duplicate_action="drop")
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = faidx.fetch('KF435150',
                         100, 150)
    assert str(result) == expect

@pytest.mark.xfail(raises=ValueError)
def test_stop(remove_index):
    fasta = Fasta('data/genes.fasta', split_char='.')

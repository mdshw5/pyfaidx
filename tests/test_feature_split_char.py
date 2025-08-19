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
    fasta = Fasta('data/genes.fasta', split_char='|', duplicate_action="drop")
    expect = ['530364724', '530364725', '530364726', '530373235', '530373237', '530384534', '530384536', '530384538', '530384540', '543583738', '543583740', '543583785', '543583786', '543583788', '543583794', '543583795', '543583796', '557361097', '557361099', '563317589', 'AB821309.1', 'KF435149.1', 'KF435150.1', 'NM_000465.3', 'NM_001282543.1', 'NM_001282545.1', 'NM_001282548.1', 'NM_001282549.1', 'NR_104212.1', 'NR_104215.1', 'NR_104216.1', 'XM_005249642.1', 'XM_005249643.1', 'XM_005249644.1', 'XM_005249645.1', 'XM_005265507.1', 'XM_005265508.1', 'XR_241079.1', 'XR_241080.1', 'XR_241081.1', 'dbj']
    result = sorted(fasta.keys())
    assert result == expect

def test_key_function_by_dictionary_get_key(remove_index):
    fasta = Fasta('data/genes.fasta', split_char='|', duplicate_action="drop")
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = fasta['KF435150.1'][100-1:150]
    assert str(result) == expect

def test_key_function_by_fetch(remove_index):
    faidx = Faidx('data/genes.fasta', split_char='|', duplicate_action="drop")
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = faidx.fetch('KF435150.1',
                         100, 150)
    assert str(result) == expect

def test_stop(remove_index):
    with pytest.raises(ValueError):
        fasta = Fasta('data/genes.fasta', split_char='|')

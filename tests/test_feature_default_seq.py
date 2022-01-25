import os
import pytest
from pyfaidx import Faidx

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_fetch_border_padded(remove_index):
    """ Fetch past the end of a gene entry """
    faidx = Faidx('data/genes.fasta', default_seq='N')
    expect = 'TCNNNNNNNNNNNNNNNNNNN'
    result = faidx.fetch('gi|557361099|gb|KF435150.1|',
                         480, 500)
    assert str(result) == expect

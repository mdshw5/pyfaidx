import os
import pytest
from pathlib import Path
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

def test_Faidx(remove_index):
    """ Ensures that Faidx can be created with a pathlib.Path as filename """
    filename = 'data/genes.fasta'
    faidx = Faidx(filename)
    faidx_w_path = Faidx(Path(filename))
    assert faidx.filename == faidx_w_path.filename

def test_Fasta(remove_index):
    """ Ensures that Fasta can be created with a pathlib.Path as filename """
    filename = 'data/genes.fasta'
    fasta = Fasta(filename)
    fasta_w_path = Fasta(Path(filename))
    assert fasta.filename == fasta_w_path.filename

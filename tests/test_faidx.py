import os
import filecmp
import pytest
from pyfaidx import BedError, FetchError
from pyfaidx.cli import main
from tempfile import NamedTemporaryFile

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file

def test_short_line_lengths(remove_index):
    with pytest.raises(BedError):
        main(['data/genes.fasta', '--bed', 'data/malformed.bed'])

def test_fetch_whole_file(remove_index):
    main(['data/genes.fasta'])

def test_split_entry(remove_index):
    main(['--split-files', 'data/genes.fasta', 'gi|557361099|gb|KF435150.1|'])
    assert os.path.exists('gi557361099gbKF435150.1.fasta')
    os.remove('gi557361099gbKF435150.1.fasta')

def test_fetch_error(remove_index):
    with pytest.raises(FetchError):
        main(['data/genes.fasta', 'gi|557361099|gb|KF435150.1|:1-1000'])
    
def test_key_warning(remove_index):
    main(['data/genes.fasta', 'foo'])
    
def test_auto_strand(remove_index):
    """ Test that --auto-strand produces the same output as --reverse --complement"""
    with NamedTemporaryFile() as auto_strand:
        with NamedTemporaryFile() as noto_strand:
            main(['--auto-strand', '-o', auto_strand.name, 'data/genes.fasta', 'gi|557361099|gb|KF435150.1|:100-1'])
            main(['--reverse', '--complement', '-o', noto_strand.name, 'data/genes.fasta', 'gi|557361099|gb|KF435150.1|:1-100'])
            print(auto_strand.read())
            print()
            print(noto_strand.read())
            assert filecmp.cmp(auto_strand.name, noto_strand.name)
    
def test_regexp(remove_index):
    main(['data/genes.fasta', '-g', 'XR'])

def test_not_regexp(remove_index):
    main(['data/genes.fasta', '-g', 'XR','-v'])

def test_not_regexp_multi(remove_index):
    main(['data/genes.fasta', '-g', 'XR', '-g', 'XM', '-v'])

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

@pytest.mark.xfail(raises=FetchError)
def test_bounds_error(remove_index):
    fasta = Fasta('data/genes.fasta', read_ahead=300, strict_bounds=True)
    result = fasta['KF435150.1'][100-1:15000].seq.lower()

@pytest.mark.xfail(raises=ValueError)
def test_buffer_value(remove_index):
    Fasta('data/genes.fasta', read_ahead=0.5)


@pytest.mark.parametrize('as_raw', [False, True])
@pytest.mark.parametrize('region, replacement, expected', [
    (1, 'T', 'ATGTACGT'),
    (slice(2, 6), 'TTAA', 'ACTTAAGT'),
])
def test_mutation_invalidates_buffer(tmp_path, as_raw, region, replacement, expected):
    filename = tmp_path / 'mutable.fa'
    filename.write_text('>chr1\nACGT\nACGT\nACGT\n')

    with Fasta(filename, mutable=True, read_ahead=4, as_raw=as_raw) as fasta:
        assert str(fasta['chr1'][:4]) == 'ACGT'

        fasta['chr1'][region] = replacement

        assert str(fasta['chr1'][:8]) == expected


@pytest.mark.parametrize('as_raw', [False, True])
def test_mutation_through_alias_invalidates_buffer(tmp_path, as_raw):
    filename = tmp_path / 'aliases.fa'
    filename.write_text('>chr1|alias\nACGT\nACGT\nACGT\n')

    with Fasta(filename, mutable=True, read_ahead=4, as_raw=as_raw,
               split_char='|') as fasta:
        assert str(fasta['chr1'][:4]) == 'ACGT'

        fasta['alias'][4:8] = 'TTAA'

        assert str(fasta['chr1'][:8]) == 'ACGTTTAA'


@pytest.mark.parametrize('as_raw', [False, True])
def test_rejected_mutation_preserves_buffered_sequence(tmp_path, as_raw):
    filename = tmp_path / 'rejected.fa'
    filename.write_text('>chr1\nACGT\nACGT\nACGT\n')
    original = filename.read_bytes()

    with Fasta(filename, mutable=True, read_ahead=4, as_raw=as_raw) as fasta:
        assert str(fasta['chr1'][:4]) == 'ACGT'

        with pytest.raises(IOError, match='same length as original'):
            fasta['chr1'][1:3] = 'N'

        assert str(fasta['chr1'][:8]) == 'ACGTACGT'
        assert filename.read_bytes() == original

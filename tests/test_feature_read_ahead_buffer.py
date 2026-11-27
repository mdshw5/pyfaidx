import os
import pytest
from pyfaidx import Faidx, Fasta, FetchError
from unittest import mock

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


@pytest.mark.parametrize('newline', [b'\n', b'\r\n'])
@pytest.mark.parametrize('as_raw', [False, True])
def test_buffer_hits_and_switches_records(tmp_path, newline, as_raw):
    filename = tmp_path / 'records.fa'
    filename.write_bytes(newline.join([
        b'>first', b'ACGT', b'TGCA', b'>second', b'GGCC', b'AATT', b''
    ]))

    with Fasta(filename, read_ahead=4, as_raw=as_raw) as fasta:
        with mock.patch.object(fasta.faidx, 'from_file', wraps=fasta.faidx.from_file) as read:
            assert str(fasta['first'][1:3]) == 'CG'
            assert str(fasta['first'][3:6]) == 'TTG'
            assert read.call_count == 1

            assert str(fasta['second'][1:3]) == 'GC'
            assert read.call_count == 2
            assert str(fasta['first'][1:3]) == 'CG'
            assert read.call_count == 3


@pytest.mark.parametrize('as_raw', [False, True])
def test_buffer_strict_bounds_near_end(tmp_path, as_raw):
    filename = tmp_path / 'record.fa'
    filename.write_text('>record\nACGT\nTGCA\n')

    with Fasta(filename, read_ahead=4, strict_bounds=True, as_raw=as_raw) as fasta:
        assert str(fasta['record'][5:8]) == 'GCA'
        with pytest.raises(FetchError):
            fasta['record'][5:9]
        assert str(fasta['record'][1:3]) == 'CG'


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

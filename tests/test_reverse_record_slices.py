import pytest

from pyfaidx import Fasta


@pytest.mark.parametrize('as_raw', [False, True])
@pytest.mark.parametrize('strict_bounds', [False, True])
@pytest.mark.parametrize('region', [
    slice(8, 1, -1), slice(8, 1, -2), slice(None, 1, -3),
    slice(7, None, -1), slice(-2, -9, -2), slice(100, -100, -1),
    slice(1, 8, -1), slice(-100, None, -1), slice(None, None, -1),
])
def test_reverse_record_slices(tmp_path, as_raw, strict_bounds, region):
    sequence = 'ACGTTAGCCTAA'
    filename = tmp_path / 'record.fa'
    filename.write_text('>record\n' + sequence + '\n')
    with Fasta(filename, as_raw=as_raw, strict_bounds=strict_bounds) as fasta:
        assert str(fasta['record'][region]) == sequence[region]


@pytest.mark.parametrize('as_raw', [False, True])
def test_reverse_empty_record(tmp_path, as_raw):
    filename = tmp_path / 'empty.fa'
    filename.write_text('>empty\n')
    with Fasta(filename, as_raw=as_raw) as fasta:
        assert str(fasta['empty'][::-1]) == ''
        assert str(fasta['empty'][5:1:-2]) == ''

import pytest
from pyfaidx import Fasta


@pytest.mark.parametrize('newline', ['\n', '\r\n'])
@pytest.mark.parametrize('start, end, expected', [
    (2, 4, ['ACnnACGTAC', 'GGTTAA']),  # within a line
    (8, 12, ['ACGTACGTnn', 'nnTTAA']),  # across a line break
    (10, 14, ['ACGTACGTAC', 'nnnnAA']),  # from the start of a line
])
def test_mutable_write(tmp_path, newline, start, end, expected):
    path = tmp_path / 'mutable.fa'
    path.write_bytes(newline.join(['>chr1', 'ACGTACGTAC', 'GGTTAA', '']).encode())
    with Fasta(str(path), mutable=True) as fa:
        fa['chr1'][start:end] = 'n' * (end - start)
    assert path.read_bytes() == newline.join(['>chr1'] + expected + ['']).encode()


@pytest.mark.parametrize('newline', ['\n', '\r\n'])
@pytest.mark.parametrize('start, end, expected', [
    (12, 20, ['ACGTACGTAC', 'GGnnnnnnnn']),  # end of the last line
    (18, 20, ['ACGTACGTAC', 'GGTTAACCnn']),
    (8, 20, ['ACGTACGTnn', 'nnnnnnnnnn']),  # across the last line break
])
def test_mutable_write_without_final_newline(tmp_path, newline, start, end, expected):
    path = tmp_path / 'mutable.fa'
    path.write_bytes(newline.join(['>chr1', 'ACGTACGTAC', 'GGTTAACCGG']).encode())
    with Fasta(str(path), mutable=True) as fa:
        fa['chr1'][start:end] = 'n' * (end - start)
    assert path.read_bytes() == newline.join(['>chr1'] + expected).encode()

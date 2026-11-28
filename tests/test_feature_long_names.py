import pytest
from pyfaidx import Fasta

Bio = pytest.importorskip("Bio.bgzf")

FASTAS = {
    'crlf': '>chr1 first record\r\nACGTACGTAC\r\nGGTTAA\r\n>chr2 second\r\nTTTTCCCCGG\r\nAA\r\n',
    'blank_line': '>chr1 first record\nACGTACGTAC\nGGTTAA\n\n>chr2 second\nTTTTCCCCGG\nAA\n',
    'comment': ';comment line\n>chr1 first record\nACGTACGTAC\nGGTTAA\n>chr2 second\nTTTTCCCCGG\nAA\n',
}


def write(tmp_path, name, text, bgzf):
    if bgzf:
        path = tmp_path / (name + '.fa.gz')
        with Bio.BgzfWriter(str(path), 'wb') as out:
            out.write(text.encode())
    else:
        path = tmp_path / (name + '.fa')
        path.write_bytes(text.encode())
    return str(path)


@pytest.mark.parametrize('bgzf', [False, True])
@pytest.mark.parametrize('name', FASTAS)
def test_long_name(tmp_path, name, bgzf):
    fa = Fasta(write(tmp_path, name, FASTAS[name], bgzf))
    assert fa['chr1'].long_name == 'chr1 first record'
    assert fa['chr2'].long_name == 'chr2 second'


@pytest.mark.parametrize('name', FASTAS)
def test_read_long_names(tmp_path, name):
    fa = Fasta(write(tmp_path, name, FASTAS[name], False), read_long_names=True)
    assert list(fa.keys()) == ['chr1 first record', 'chr2 second']

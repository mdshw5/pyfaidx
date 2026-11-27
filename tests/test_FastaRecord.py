import os
import sys
import pytest
from pyfaidx import Fasta
from tempfile import NamedTemporaryFile
from difflib import Differ

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file
    
def test_sequence_uppercase():
    """Test that the sequence is always returned in
    uppercase, even if it is in lowercase in the
    reference genome.
    """
    filename = "data/genes.fasta.lower"
    reference_upper = Fasta(filename, sequence_always_upper=True)
    reference_normal = Fasta(filename)
    os.remove('data/genes.fasta.lower.fai')
    assert reference_upper['KF435150.1'][
    1:100].seq == reference_normal['KF435150.1'][
            1:100].seq.upper()

def test_long_names():
    """ Test that deflines extracted using FastaRecord.long_name are
    identical to deflines in the actual file.
    """
    deflines = []
    with open('data/genes.fasta') as fasta_file:
        for line in fasta_file:
            if line[0] == '>':
                deflines.append(line[1:-1])
    fasta = Fasta('data/genes.fasta')
    long_names = []
    for record in fasta:
        long_names.append(record.long_name)
    print(tuple(zip(deflines, long_names)))
    assert deflines == long_names

def write_uniform_fasta(fasta_file, fasta_uniform_len):
    deflines = []
    line_len = None
    newline = None
    for line in fasta_file:
        if line.startswith(b'>'):
            deflines.append(line[1:].rstrip(b'\r\n').decode('ascii'))
            fasta_uniform_len.write(line)
        elif line_len is None:
            line_len = len(line)
            newline = b'\r\n' if line.endswith(b'\r\n') else b'\n'
            fasta_uniform_len.write(line)
        elif line_len > len(line):
            fasta_uniform_len.write(line.rstrip(b'\r\n') + b'N' *
                                    (line_len - len(line)) + newline)
        else:
            fasta_uniform_len.write(line)
    return deflines

def test_issue_62():
    """ Check for pathogenic FastaRecord.long_name behavior in mdshw5/pyfaidx#62 """
    with open('data/genes.fasta', 'rb') as fasta_file:
        with open('data/issue_62.fa', 'wb') as fasta_uniform_len:
            deflines = write_uniform_fasta(fasta_file, fasta_uniform_len)
    fasta = Fasta('data/issue_62.fa', as_raw=True)
    long_names = []
    for record in fasta:
        long_names.append(record.long_name)
    try:
        os.remove('data/issue_62.fa')
        os.remove('data/issue_62.fa.fai')
    except EnvironmentError:
        pass
    sys.stdout.writelines(tuple(Differ().compare(deflines, long_names)))
    assert deflines == long_names

@pytest.mark.parametrize('newline', [b'\n', b'\r\n'])
def test_issue_62_uniform_line_endings(tmp_path, newline):
    from io import BytesIO
    source = BytesIO(newline.join([b'>first description', b'ACGT', b'AC',
                                   b'>second description', b'TGCA', b'TG', b'']))
    filename = tmp_path / 'uniform.fa'
    with filename.open('wb') as output:
        expected = write_uniform_fasta(source, output)
    assert filename.read_bytes() == newline.join([
        b'>first description', b'ACGT', b'ACNN',
        b'>second description', b'TGCA', b'TGNN', b''])
    with Fasta(str(filename), as_raw=True) as fasta:
        assert [record.long_name for record in fasta] == expected

def test_unpadded_length():
    filename = "data/padded.fasta"
    with open(filename, 'w') as padded:
        padded.write(">test_padded\n")
        for n in range(10):
            padded.write("N" * 80)
            padded.write("\n")
        padded.write("N" * 30)
        padded.write("A" * 20)
        padded.write("N" * 30)
        padded.write("\n")
        for n in range(10):
            padded.write("N" * 80)
            padded.write("\n")

    fasta = Fasta(filename)
    expect = 20
    result = fasta["test_padded"].unpadded_len
    print(expect, result)
    assert expect == result
    fasta.close()
    os.remove('data/padded.fasta')
    os.remove('data/padded.fasta.fai')

@pytest.mark.parametrize('as_raw', [False, True])
def test_unpadded_length_all_n_and_empty(tmp_path, as_raw):
    filename = tmp_path / 'all_n_and_empty.fasta'
    filename.write_text(
        '>all_n\n' + ('N' * 80 + '\n') * 2 +
        '>all_n_lower\n' + ('n' * 50 + '\n') * 2 +
        '>empty\n' +
        '>single_base\n' + 'N' * 10 + 'A' + 'N' * 10 + '\n' +
        '>padded\n' + 'N' * 5 + 'ACGT' + 'N' * 5 + '\n'
    )
    with Fasta(str(filename), as_raw=as_raw) as fasta:
        assert fasta['all_n'].unpadded_len == 0
        assert fasta['all_n_lower'].unpadded_len == 0
        assert fasta['empty'].unpadded_len == 0
        assert fasta['single_base'].unpadded_len == 1
        assert fasta['padded'].unpadded_len == 4

def test_numpy_array(remove_index):
    """ Test the __array_interface__ """
    import numpy
    filename = "data/genes.fasta.lower"
    reference = Fasta(filename)
    np_array = numpy.asarray(reference[0])
    assert isinstance(np_array, numpy.ndarray)

@pytest.fixture
def remove_index_mutable():
    with open('data/genes_mutable.fasta', 'wb') as mutable:
        mutable.write(open('data/genes.fasta', 'rb').read())
    mutable_fasta = Fasta('data/genes_mutable.fasta', mutable=True)
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file
    try:
        os.remove('data/genes_mutable.fasta')
    except EnvironmentError:
        pass  # some tests may delete this file
    try:
        os.remove('data/genes_mutable.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file


@pytest.mark.usefixtures('remove_index_mutable')
def test_mutate_fasta_to_same():
    mutable = Fasta('data/genes_mutable.fasta', mutable=True)
    fasta = Fasta('data/genes.fasta', mutable=False)
    chunk = fasta['KF435150.1'][0:100]
    mutable['KF435150.1'][0:100] = chunk.seq
    assert str(fasta['KF435150.1']) == str(mutable['KF435150.1'])

@pytest.mark.usefixtures('remove_index_mutable')
def test_mutate_fasta_to_N():
    mutable = Fasta('data/genes_mutable.fasta', mutable=True)
    chunk = 100 * 'N'
    mutable['KF435150.1'][0:100] = chunk
    assert mutable['KF435150.1'][0:100].seq == chunk

def test_mutate_crlf_fasta_preserves_line_endings(tmp_path):
    filename = tmp_path / 'crlf.fasta'
    filename.write_bytes(b'>chr1\r\nACGT\r\nACGT\r\n')

    with Fasta(str(filename), mutable=True) as fasta:
        fasta['chr1'][1:7] = 'CCCCCC'

    assert filename.read_bytes() == b'>chr1\r\nACCC\r\nCCCT\r\n'

@pytest.mark.usefixtures('remove_index_mutable')
def test_mutate_single_position():
    mutable = Fasta('data/genes_mutable.fasta', mutable=True)
    chunk = 'N'
    mutable['KF435150.1'][0] = chunk
    assert mutable['KF435150.1'][0].seq == chunk

@pytest.mark.xfail(raises=TypeError)
@pytest.mark.usefixtures('remove_index_mutable')
def test_mutate_immutable_fasta():
    mutable = Fasta('data/genes_mutable.fasta', mutable=False)
    chunk = 100 * 'N'
    mutable['KF435150.1'][0:100] = chunk

@pytest.mark.xfail(raises=IOError)
@pytest.mark.usefixtures('remove_index_mutable')
def test_mutate_too_long():
    mutable = Fasta('data/genes_mutable.fasta', mutable=True)
    chunk = 101 * 'N'
    mutable['KF435150.1'][0:100] = chunk

def test_fast_sum_FastaRecord_len():
    """ Test that the sum of FastaRecord lengths is correct """
    fasta = Fasta('data/genes.fasta')
    total_length = sum(len(record) for record in fasta)
    assert total_length == len(fasta) 

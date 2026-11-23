import os
import filecmp
import pytest
from pyfaidx import BedError, FetchError, FastaNotFoundError, Fasta, Faidx
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

def test_short_line_lengths():
    with pytest.raises(BedError):
        main(['data/genes.fasta', '--bed', 'data/malformed.bed'])

def test_fetch_whole_file():
    main(['data/genes.fasta'])

def test_split_entry():
    main(['--split-files', 'data/genes.fasta', 'KF435150.1'])
    assert os.path.exists('KF435150.1.fasta')
    os.remove('KF435150.1.fasta')

@pytest.mark.xfail(raises=FetchError)
def test_fetch_error():
    main(['data/genes.fasta', 'KF435150.1:1-1000'])
    
def test_key_warning():
    main(['data/genes.fasta', 'foo'])
    
def test_auto_strand(tmp_path):
    """ Test that --auto-strand produces the same output as --reverse --complement"""
    fasta = tmp_path / 'genes.fasta'
    fasta.write_text('>KF435150.1\n' + ('AACCGTTA' * 50) + '\n')
    auto_strand = tmp_path / 'auto-strand.fa'
    noto_strand = tmp_path / 'reverse-complement.fa'
    main(['--auto-strand', '--no-names', '-o', str(auto_strand), str(fasta), 'KF435150.1:100-1'])
    main(['--reverse', '--complement', '--no-names', '-o', str(noto_strand), str(fasta), 'KF435150.1:1-100'])
    assert filecmp.cmp(auto_strand, noto_strand)


@pytest.mark.parametrize('regions, expected', [
    ([], 'ACGTAC\nGATT\n'),
    (['test'], 'ACGTAC\n'),
    (['test:2-5'], 'CGTA\n'),
    (['test:5-2'], 'TACG\n'),
])
def test_auto_strand_regions(tmp_path, capsys, regions, expected):
    filename = tmp_path / 'reference.fa'
    filename.write_bytes(b'>test\nACGTAC\n>other\nGATT\n')
    main(['--auto-strand', '--no-names', str(filename)] + regions)
    assert capsys.readouterr().out == expected
    
def test_regexp():
    main(['data/genes.fasta', '-g', 'XR'])

def test_not_regexp():
    main(['data/genes.fasta', '-g', 'XR','-v'])

def test_not_regexp_multi():
    main(['data/genes.fasta', '-g', 'XR', '-g', 'XM', '-v'])

@pytest.mark.xfail(raises=SystemExit)
def test_faidx_cli_invalid_args():
    """Test that the CLI raises an error for invalid arguments."""
    main(['--invalid-arg'])

@pytest.mark.xfail(raises=FastaNotFoundError)
def test_faidx_cli_missing_file():
    """Test that the CLI raises an error for a missing file."""
    main(['nonexistent_file.fasta'])

def test_faidx_cli_transform_nucleotides():
    """Test that the CLI can transform nucleotides."""
    with NamedTemporaryFile(delete=False) as temp_fasta:
        with NamedTemporaryFile(delete=False) as transform_table:
            temp_fasta.write(b'>test\nATCG\n')
            temp_fasta.close()
            result = main(['--out', transform_table.name, '--transform', 'nucleotide', temp_fasta.name ])
            assert result is None  # CLI should not return anything
            with open(transform_table.name, 'r') as f:
                content = f.read()
            # Check that the output matches the expected format
            assert content == 'name\tstart\tend\tA\tT\tC\tG\tN\tothers\n' +\
                            'test\t1\t4\t1\t1\t1\t1\t0\t0\n'
    # Clean up temporary files
    os.remove(temp_fasta.name)
    os.remove(transform_table.name)

def test_faidx_cli_transform_transposed():
    """Test that the CLI can transform sequences to transposed format."""
    with NamedTemporaryFile(delete=False) as temp_fasta:
        with NamedTemporaryFile(delete=False) as transform_table:
            temp_fasta.write(b'>test\nATCG\n')
            temp_fasta.close()
            result = main(['--out', transform_table.name, '--transform', 'transposed', temp_fasta.name])
            assert result is None  # CLI should not return anything
            with open(transform_table.name, 'r') as f:
                content = f.read()
            # Check that the output matches the expected format
            assert content == 'test\t1\t4\tATCG\n'
    # Clean up temporary files
    os.remove(temp_fasta.name)
    os.remove(transform_table.name)

@pytest.mark.parametrize('transform, expected', [
    ('transposed', 'test\t4\t1\tGCAA\n'),
    ('nucleotide', 'name\tstart\tend\tA\tT\tC\tG\tN\tothers\n'
                   'test\t4\t1\t2\t0\t1\t1\t0\t0\n'),
])
def test_faidx_cli_transform_reverse(tmp_path, capsys, transform, expected):
    fasta = tmp_path / 'reference.fa'
    fasta.write_bytes(b'>test\nAACG\n')
    main(['--reverse', '--transform', transform, str(fasta)])
    assert capsys.readouterr().out == expected


def test_faidx_cli_transform_bed():
    """Test that the CLI can transform sequences to BED format."""
    with NamedTemporaryFile(delete=False) as temp_fasta:
        with NamedTemporaryFile(delete=False) as transform_bed:
            temp_fasta.write(b'>test\nATCG\n')
            temp_fasta.close()
            result = main(['--out', transform_bed.name, '--transform', 'bed', temp_fasta.name])
            assert result is None  # CLI should not return anything
            with open(transform_bed.name, 'r') as f:
                content = f.read()
            # Check that the output matches the expected format
            assert content == 'test\t0\t4\n'
    # Clean up temporary files
    os.remove(temp_fasta.name)
    os.remove(transform_bed.name)

@pytest.mark.parametrize('compressed', [False, True])
def test_faidx_cli_transform_chromsizes(tmp_path, capsys, monkeypatch, compressed):
    """Whole-record lengths should not require fetching chromosome sequences."""
    filename = tmp_path / ('chromosomes.fa.gz' if compressed else 'chromosomes.fa')
    content = b'>chr1 description\nATCG\nATCG\n>chr2\nATCGA\n'
    if compressed:
        from Bio.bgzf import BgzfWriter
        with BgzfWriter(str(filename), 'wb') as fasta:
            fasta.write(content)
    else:
        filename.write_bytes(content)
    with Fasta(str(filename)):
        pass  # Build the indexes before checking indexed retrieval.

    def unexpected_fetch(*args, **kwargs):
        pytest.fail('chromsizes fetched sequence instead of using the index')

    monkeypatch.setattr(Faidx, 'fetch', unexpected_fetch)
    main(['--transform', 'chromsizes', str(filename)])
    assert capsys.readouterr().out == 'chr1\t8\nchr2\t5\n'

def test_faidx_cli_size_range():
    """Test that the CLI can handle size range arguments."""
    with NamedTemporaryFile(delete=False) as temp_fasta:
        with NamedTemporaryFile(delete=False) as output_file:
            temp_fasta.write(b'>test1\nATCG\n' + b'>test2\nAT\n' + b'>test3\nATCGTAGC\n')
            temp_fasta.close()
            result = main(['--size-range', '3,5', '--out', output_file.name, temp_fasta.name])
            assert result is None  # CLI should not return anything
            # Check that the output contains only sequences within the size range
            with open(output_file.name, 'r') as f:
                content = f.read()
            assert '>test1' in content
            assert '>test2' not in content
            assert '>test3' not in content
    # Clean up temporary files
    os.remove(temp_fasta.name)
    os.remove(output_file.name)

def test_default_seq_length():
    """Test that passing a default sequence > 1 raises an error."""
    with pytest.raises(SystemExit):
        main(['data/genes.fasta', '--default-seq', 'NN'])

def test_faidx_exit_on_missing_file(remove_index):
    """Test that Faidx.__exit__ does not raise AttributeError if file is missing. (#229)"""
    class FaidxWrapper:
        def __init__(self):
            try:
                self.obj = Faidx('nonexistent_file.fasta')
            except Exception:
                pass
        def __del__(self):
            # Should not raise AttributeError
            if hasattr(self, 'obj'):
                del self.obj


    wrapper = FaidxWrapper()
    # Deleting wrapper should not raise AttributeError
    try:
        del wrapper
    except AttributeError as e:
        pytest.fail(f"Unexpected AttributeError on __del__: {e}")

def test_faidx_sequence_masking_default():
    import time
    with NamedTemporaryFile(delete=False) as temp_fasta:
        temp_fasta.write(b'>test\nATCG\n')
        temp_fasta_name = temp_fasta.name
    result = main(['--mask-with-default-seq', '--default-seq', 'N', temp_fasta_name, 'test:1-2'])
    assert result is None  # CLI should not return anything
    with open(temp_fasta_name, 'r') as f:
        content = f.read()
    # Check that the output matches the expected format
    assert content == '>test\nNNCG\n'
    # Clean up temporary files
    os.remove(temp_fasta_name)

def test_faidx_sequence_masking_case():
    with NamedTemporaryFile(delete=False) as temp_fasta:
        temp_fasta.write(b'>test\nATCG\n')
        temp_fasta_name = temp_fasta.name
    result = main(['--mask-by-case', temp_fasta_name, 'test:1-2'])
    assert result is None  # CLI should not return anything
    with open(temp_fasta_name, 'r') as f:
        content = f.read()
    # Check that the output matches the expected format
    assert content == '>test\natCG\n'
    # Clean up temporary files
    os.remove(temp_fasta_name)


@pytest.mark.parametrize('action, selected', [
    ('first', 0), ('last', 2), ('longest', 1), ('shortest', 2)])
@pytest.mark.parametrize('mode', ['extract', 'transform', 'soft-mask', 'hard-mask'])
def test_faidx_cli_duplicates_action(tmp_path, capsys, action, selected, mode):
    sequences = ['AAAA', 'CCCCCC', 'TT']
    filename = tmp_path / 'duplicates.fa'
    filename.write_bytes(''.join('>chr1\n' + seq + '\n' for seq in sequences).encode())
    args = [str(filename)]
    if mode == 'extract':
        args += ['--no-names']
    elif mode == 'transform':
        args += ['--transform', 'transposed']
    elif mode == 'soft-mask':
        args += ['chr1:1-2', '--mask-by-case']
    else:
        args += ['chr1:1-2', '--mask-with-default-seq', '--default-seq', 'N']
    args += ['--duplicates-action', action]

    main(args)
    output = capsys.readouterr().out
    sequence = sequences[selected]
    if mode == 'extract':
        assert output == sequence + '\n'
    elif mode == 'transform':
        assert output == 'chr1\t1\t{}\t{}\n'.format(len(sequence), sequence)
    else:
        assert output == ''
        prefix = sequence[:2].lower() if mode == 'soft-mask' else 'NN'
        sequences[selected] = prefix + sequence[2:]
    assert filename.read_text() == ''.join('>chr1\n' + seq + '\n' for seq in sequences)


@pytest.mark.parametrize('action', [[], ['--duplicates-action', 'stop']])
def test_faidx_cli_duplicate_names_stop(tmp_path, action):
    filename = tmp_path / 'duplicates.fa'
    filename.write_bytes(b'>chr1\nAAAA\n>chr1\nTT\n')
    with pytest.raises(ValueError, match='Duplicate key "chr1"'):
        main([str(filename)] + action)

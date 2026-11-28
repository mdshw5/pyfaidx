import pytest

from pyfaidx import BedError, bed_split
from pyfaidx.cli import main


@pytest.mark.parametrize('entry, expected', [
    ('chr1\t000\t10', ('chr1', 0, 10)),
    ('chr1 10 10', ('chr1', 10, 10)),
    ('chr1\t0\t0\tfeature', ('chr1', 0, 0)),
])
def test_bed_split_valid_entries(entry, expected):
    assert bed_split(entry) == expected


@pytest.mark.parametrize('entry', [
    'chr1\t1.5\t2',
    'chr1\tstart\t2',
    'chr1\t-1\t2',
    'chr1\t1\t-2',
    'chr1\t3\t2',
    'chr1\t0',
    '\t0\t2',
    'chr1&\t0\t2',
    '"chr1"\t0\t2',
    'chr1\t' + '9' * 5000 + '\t1',
])
def test_bed_split_rejects_malformed_entries(entry):
    with pytest.raises(BedError):
        bed_split(entry)


@pytest.mark.parametrize('non_data', ['# comment\n', '\n', ' \t \n', '  # comment\n'])
@pytest.mark.parametrize('mode', ['extract', 'transform', 'soft-mask', 'hard-mask'])
def test_bed_non_data_lines(tmp_path, capsys, non_data, mode):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\nACGTACGT\n')
    bed = tmp_path / 'regions.bed'
    bed.write_text(non_data + 'chr1\t0\t2\n' + non_data + 'chr1\t4\t6\n' + non_data)
    args = [str(fasta), '--bed', str(bed)]
    if mode == 'extract':
        args += ['--no-names']
    elif mode == 'transform':
        args += ['--transform', 'bed']
    elif mode == 'soft-mask':
        args += ['--mask-by-case']
    else:
        args += ['--mask-with-default-seq', '--default-seq', 'N']
    main(args)
    output = capsys.readouterr().out
    output = output.replace('\r\n', '\n')
    if mode == 'extract':
        assert output == 'AC\nAC\n'
    elif mode == 'transform':
        assert output == 'chr1\t0\t2\nchr1\t4\t6\n'
    else:
        assert output == ''
    expected = {'soft-mask': 'acGTacGT', 'hard-mask': 'NNGTNNGT'}.get(mode, 'ACGTACGT')
    assert fasta.read_text() == '>chr1\n' + expected + '\n'


@pytest.mark.parametrize('mask', [[], ['--mask-by-case'], ['--mask-with-default-seq', '--default-seq', 'N']])
def test_bed_with_no_data_leaves_reference_unchanged(tmp_path, capsys, mask):
    fasta = tmp_path / 'reference.fa'
    original = '>chr1\nACGT\n'
    fasta.write_text(original)
    bed = tmp_path / 'empty.bed'
    bed.write_text('\n# no regions\n\t\n')
    main([str(fasta), '--bed', str(bed)] + mask)
    assert capsys.readouterr().out == ''
    assert fasta.read_text() == original


def test_nonempty_malformed_bed_still_fails(tmp_path):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\nACGT\n')
    bed = tmp_path / 'bad.bed'
    bed.write_text('chr1\t0\n')
    with pytest.raises(BedError):
        main([str(fasta), '--bed', str(bed)])


def test_zero_length_bed_interval_transforms_correctly(tmp_path, capsys):
    fasta = tmp_path / 'reference.fa'
    fasta.write_text('>chr1\nACGT\n')
    bed = tmp_path / 'zero-length.bed'
    bed.write_text('chr1\t0\t0\n')

    main([str(fasta), '--bed', str(bed), '--transform', 'bed'])

    assert capsys.readouterr().out == 'chr1\t0\t0\n'


def test_zero_length_bed_interval_does_not_mask_reference(tmp_path):
    fasta = tmp_path / 'reference.fa'
    original = '>chr1\nACGT\n'
    fasta.write_text(original)
    bed = tmp_path / 'zero-length.bed'
    bed.write_text('chr1\t0\t0\n')

    main([str(fasta), '--bed', str(bed), '--mask-with-default-seq',
          '--default-seq', 'N'])

    assert fasta.read_text() == original

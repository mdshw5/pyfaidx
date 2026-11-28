import os
import pytest
from pyfaidx import Fasta, FastaIndexingError


def write(tmp_path, text, name='comments.fasta'):
    path = tmp_path / name
    path.write_bytes(text.encode())
    return str(path)


def test_leading_comment_lines_are_skipped(tmp_path):
    fasta = write(tmp_path, ";comment line\n;~key: value\n;\n>chr1 first\nACGTACGT\nACGT\n>chr2\nTTTTGGGG\n")
    fa = Fasta(fasta)
    assert list(fa.keys()) == ['chr1', 'chr2']
    assert str(fa['chr1'][:]) == 'ACGTACGTACGT'
    assert str(fa['chr1'][2:10]) == 'GTACGTAC'
    assert str(fa['chr2'][:]) == 'TTTTGGGG'
    with open(fasta + '.fai') as fai:
        assert fai.read() == "chr1\t12\t41\t8\t9\nchr2\t8\t61\t8\t9\n"


def test_leading_comment_lines_of_different_lengths(tmp_path):
    fasta = write(tmp_path, ";a longer comment line\n;short\n>chr1\nACGT\n")
    fa = Fasta(fasta)
    assert list(fa.keys()) == ['chr1']
    assert str(fa['chr1'][:]) == 'ACGT'


def test_leading_comment_and_blank_lines_crlf(tmp_path):
    fasta = write(tmp_path, ";comment\r\n\r\n>chr1\r\nACGT\r\nAC\r\n")
    fa = Fasta(fasta)
    assert list(fa.keys()) == ['chr1']
    assert str(fa['chr1'][:]) == 'ACGTAC'


def test_text_before_first_record_is_an_error(tmp_path):
    fasta = write(tmp_path, "not a comment\n>chr1\nACGT\n")
    with pytest.raises(FastaIndexingError, match="before the first sequence definition line"):
        Fasta(fasta)


def test_comments_only_is_an_error(tmp_path):
    fasta = write(tmp_path, ";only a comment\n")
    with pytest.raises(FastaIndexingError):
        Fasta(fasta)

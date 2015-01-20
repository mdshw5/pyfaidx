import os
from pyfaidx import FastaIndexingError, BedError, FetchError
from pyfaidx.cli import main
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)


class TestCLI:

    @raises(FastaIndexingError)
    def test_short_line_lengths(self):
        main(['data/short_line.fasta'])

    @raises(BedError)
    def test_short_line_lengths(self):
        main(['data/genes.fasta', '--bed', 'data/malformed.bed'])

    def test_fetch_whole_file(self):
        main(['data/genes.fasta'])

    def test_split_entry(self):
        main(['--split-files', 'data/genes.fasta', 'KF435150.1'])
        assert os.path.exists('KF435150.1.fasta')
        os.remove('KF435150.1.fasta')

    @raises(FetchError)
    def test_fetch_error(self):
        main(['data/genes.fasta', 'KF435150.1:1-1000'])

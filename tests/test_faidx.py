import os
from pyfaidx import FastaIndexingError, BedError
from pyfaidx.cli import main
from nose.tools import raises

path = os.path.dirname(__file__)
os.chdir(path)


class TestCLI:

    def teardownclass(self):
        os.remove('data/short_line.fasta.fai')

    @raises(FastaIndexingError)
    def test_short_line_lengths(self):
        main(['data/short_line.fasta'])

    @raises(BedError)
    def test_short_line_lengths(self):
        main(['data/genes.fasta', '--bed', 'data/malformed.bed'])

import os
from unittest import TestCase
from pathlib import Path

from pyfaidx import Faidx, Fasta


path = os.path.dirname(__file__)
os.chdir(path)


class TestAcceptPath(TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        try:
            os.remove('data/genes.fasta.fai')
        except EnvironmentError:
            pass  # some tests may delete this file

    def test_Faidx(self):
        """ Ensures that Faidx can be created with a pathlib.Path as filename """
        filename = 'data/genes.fasta'
        faidx = Faidx(filename)
        faidx_w_path = Faidx(Path(filename))
        assert faidx.filename == faidx_w_path.filename

    def test_Fasta(self):
        """ Ensures that Fasta can be created with a pathlib.Path as filename """
        filename = 'data/genes.fasta'
        fasta = Fasta(filename)
        fasta_w_path = Fasta(Path(filename))
        assert fasta.filename == fasta_w_path.filename









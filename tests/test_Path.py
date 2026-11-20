import os
from pathlib import Path
from pyfaidx import Faidx, Fasta

def test_Faidx(tmp_path):
    """ Ensures that Faidx can be created with a pathlib.Path as filename """
    filename = tmp_path / 'genes.fasta'
    filename.write_text('>chr1\nACGT\n')
    with Faidx(str(filename)) as faidx, Faidx(Path(filename)) as faidx_w_path:
        assert os.path.samefile(faidx.filename, faidx_w_path.filename)

def test_Fasta(tmp_path):
    """ Ensures that Fasta can be created with a pathlib.Path as filename """
    filename = tmp_path / 'genes.fasta'
    filename.write_text('>chr1\nACGT\n')
    with Fasta(str(filename)) as fasta, Fasta(Path(filename)) as fasta_w_path:
        assert os.path.samefile(fasta.filename, fasta_w_path.filename)

import os

import pytest

from pyfaidx import Fasta

try:
    import fsspec
    from fsspec.core import OpenFile
except ImportError:
    pytestmark = pytest.mark.skip


@pytest.fixture(scope="session")
def openfile_genes_fasta():
    testdir = os.path.dirname(__file__)
    genes_fasta = os.path.join(testdir, 'data', 'genes.fasta')

    fs = fsspec.get_filesystem_class("memory")()
    with fs.open("genes.fasta", mode="wb") as f:
        with open(genes_fasta, mode="rb") as g:
            f.write(g.read())

    yield OpenFile(fs, 'genes.fasta', mode='rb')



def test_fetch_whole_file(openfile_genes_fasta):
    _ = Fasta(openfile_genes_fasta)

import os

import pytest

from pyfaidx import Fasta

try:
    import fsspec
    from fsspec.core import OpenFile
except ImportError:
    pytestmark = pytest.mark.skip


@pytest.fixture(scope="function")
def openfile_genes_fasta():
    testdir = os.path.dirname(__file__)
    genes_fasta = os.path.join(testdir, 'data', 'genes.fasta')

    fs = fsspec.filesystem("memory")
    with fs.open('genes.fasta', mode='wb') as f:
        with open(genes_fasta, mode="rb") as g:
            f.write(g.read())

    try:
        yield fsspec.open('memory://genes.fasta', mode='rb')
    finally:
        fs.rm("/**", recursive=True)
        assert not fs.ls("/")


def test_fsspec_fetch_whole_file(openfile_genes_fasta):
    _ = Fasta(openfile_genes_fasta)


def test_fsspec_default_index(openfile_genes_fasta):
    _ = Fasta(openfile_genes_fasta)

    fs = openfile_genes_fasta.fs
    assert fs.isfile(openfile_genes_fasta.path + ".fai")
    assert fs.size(openfile_genes_fasta.path + ".fai") > 0


def test_fsspec_local_index(openfile_genes_fasta, tmp_path):
    index = tmp_path.joinpath("my_local_index.fai")
    _ = Fasta(openfile_genes_fasta, indexname=index)
    assert index.is_file()
    assert index.stat().st_size > 0


def test_fsspec_remote_index(openfile_genes_fasta):
    f_fai = fsspec.open("memory://some_other_index.fai")
    _ = Fasta(openfile_genes_fasta, indexname=f_fai)

    fs = f_fai.fs
    assert fs.isfile("some_other_index.fai")
    assert fs.size("some_other_index.fai") > 0

import builtins
import os
import pytest
from os.path import getmtime
from pyfaidx import Faidx, FastaIndexingError, IndexNotFoundError, FastaNotFoundError
from tempfile import NamedTemporaryFile, mkdtemp
import time
import platform
import shutil

from unittest import mock

path = os.path.dirname(__file__)
os.chdir(path)

@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.fai')
    except EnvironmentError:
        pass  # some tests may delete this file
    
def test_version_issue_206():
    import pyfaidx
    assert isinstance(pyfaidx.__version__, str)    

def test_build(remove_index):
    expect_index = ("gi|563317589|dbj|AB821309.1|	3510	114	70	71\n"
                    "gi|557361099|gb|KF435150.1|	481	3789	70	71\n"
                    "gi|557361097|gb|KF435149.1|	642	4368	70	71\n"
                    "gi|543583796|ref|NR_104216.1|	4573	5141	70	71\n"
                    "gi|543583795|ref|NR_104215.1|	5317	9901	70	71\n"
                    "gi|543583794|ref|NR_104212.1|	5374	15415	70	71\n"
                    "gi|543583788|ref|NM_001282545.1|	4170	20980	70	71\n"
                    "gi|543583786|ref|NM_001282543.1|	5466	25324	70	71\n"
                    "gi|543583785|ref|NM_000465.3|	5523	30980	70	71\n"
                    "gi|543583740|ref|NM_001282549.1|	3984	36696	70	71\n"
                    "gi|543583738|ref|NM_001282548.1|	4113	40851	70	71\n"
                    "gi|530384540|ref|XM_005249645.1|	2752	45151	70	71\n"
                    "gi|530384538|ref|XM_005249644.1|	3004	48071	70	71\n"
                    "gi|530384536|ref|XM_005249643.1|	3109	51246	70	71\n"
                    "gi|530384534|ref|XM_005249642.1|	3097	54528	70	71\n"
                    "gi|530373237|ref|XM_005265508.1|	2794	57830	70	71\n"
                    "gi|530373235|ref|XM_005265507.1|	2848	60824	70	71\n"
                    "gi|530364726|ref|XR_241081.1|	1009	63849	70	71\n"
                    "gi|530364725|ref|XR_241080.1|	4884	65009	70	71\n"
                    "gi|530364724|ref|XR_241079.1|	2819	70099	70	71\n")
    index_file = Faidx('data/genes.fasta').indexname
    result_index = open(index_file).read()
    assert result_index == expect_index

def test_build_issue_141(remove_index):
    expect_index = ("gi|563317589|dbj|AB821309.1|	3510	115	70	72\n"
                    "gi|557361099|gb|KF435150.1|	481	3842	70	72\n"
                    "gi|557361097|gb|KF435149.1|	642	4429	70	72\n"
                    "gi|543583796|ref|NR_104216.1|	4573	5213	70	72\n"
                    "gi|543583795|ref|NR_104215.1|	5317	10040	70	72\n"
                    "gi|543583794|ref|NR_104212.1|	5374	15631	70	72\n"
                    "gi|543583788|ref|NM_001282545.1|	4170	21274	70	72\n"
                    "gi|543583786|ref|NM_001282543.1|	5466	25679	70	72\n"
                    "gi|543583785|ref|NM_000465.3|	5523	31415	70	72\n"
                    "gi|543583740|ref|NM_001282549.1|	3984	37211	70	72\n"
                    "gi|543583738|ref|NM_001282548.1|	4113	41424	70	72\n"
                    "gi|530384540|ref|XM_005249645.1|	2752	45784	70	72\n"
                    "gi|530384538|ref|XM_005249644.1|	3004	48745	70	72\n"
                    "gi|530384536|ref|XM_005249643.1|	3109	51964	70	72\n"
                    "gi|530384534|ref|XM_005249642.1|	3097	55292	70	72\n"
                    "gi|530373237|ref|XM_005265508.1|	2794	58640	70	72\n"
                    "gi|530373235|ref|XM_005265507.1|	2848	61675	70	72\n"
                    "gi|530364726|ref|XR_241081.1|	1009	64742	70	72\n"
                    "gi|530364725|ref|XR_241080.1|	4884	65918	70	72\n"
                    "gi|530364724|ref|XR_241079.1|	2819	71079	70	72\n")
    index_file = Faidx('data/issue_141.fasta').indexname
    result_index = open(index_file).read()
    os.remove('data/issue_141.fasta.fai')
    print(result_index)
    assert result_index == expect_index

def test_build_issue_111(remove_index):
    expect_index = ("gi|563317589|dbj|AB821309	3510	114	70	71\n"
                    "gi|557361099|gb|KF435150	481	3789	70	71\n"
                    "gi|557361097|gb|KF435149	642	4368	70	71\n"
                    "gi|543583796|ref|NR_104216	4573	5141	70	71\n"
                    "gi|543583795|ref|NR_104215	5317	9901	70	71\n"
                    "gi|543583794|ref|NR_104212	5374	15415	70	71\n"
                    "gi|543583788|ref|NM_001282545	4170	20980	70	71\n"
                    "gi|543583786|ref|NM_001282543	5466	25324	70	71\n"
                    "gi|543583785|ref|NM_000465	5523	30980	70	71\n"
                    "gi|543583740|ref|NM_001282549	3984	36696	70	71\n"
                    "gi|543583738|ref|NM_001282548	4113	40851	70	71\n"
                    "gi|530384540|ref|XM_005249645	2752	45151	70	71\n"
                    "gi|530384538|ref|XM_005249644	3004	48071	70	71\n"
                    "gi|530384536|ref|XM_005249643	3109	51246	70	71\n"
                    "gi|530384534|ref|XM_005249642	3097	54528	70	71\n"
                    "gi|530373237|ref|XM_005265508	2794	57830	70	71\n"
                    "gi|530373235|ref|XM_005265507	2848	60824	70	71\n"
                    "gi|530364726|ref|XR_241081	1009	63849	70	71\n"
                    "gi|530364725|ref|XR_241080	4884	65009	70	71\n"
                    "gi|530364724|ref|XR_241079	2819	70099	70	71\n")
    index = Faidx(
        'data/genes.fasta',
        read_long_names=True,
        key_function=lambda x: x.split('.')[0])
    result_index = ''.join(index._index_as_string())
    assert result_index == expect_index

def test_order(remove_index):
    order = ("gi|563317589|dbj|AB821309.1|", "gi|557361099|gb|KF435150.1|",
             "gi|557361097|gb|KF435149.1|",
             "gi|543583796|ref|NR_104216.1|",
             "gi|543583795|ref|NR_104215.1|",
             "gi|543583794|ref|NR_104212.1|",
             "gi|543583788|ref|NM_001282545.1|",
             "gi|543583786|ref|NM_001282543.1|",
             "gi|543583785|ref|NM_000465.3|",
             "gi|543583740|ref|NM_001282549.1|",
             "gi|543583738|ref|NM_001282548.1|",
             "gi|530384540|ref|XM_005249645.1|",
             "gi|530384538|ref|XM_005249644.1|",
             "gi|530384536|ref|XM_005249643.1|",
             "gi|530384534|ref|XM_005249642.1|",
             "gi|530373237|ref|XM_005265508.1|",
             "gi|530373235|ref|XM_005265507.1|",
             "gi|530364726|ref|XR_241081.1|",
             "gi|530364725|ref|XR_241080.1|",
             "gi|530364724|ref|XR_241079.1|")
    result = tuple(Faidx('data/genes.fasta').index.keys())
    assert result == order

def test_valgrind_short_lines(remove_index):
    """ Makes all full-length lines short and checks that error is raised
    in all appropriate circumstances.
    """
    # http://stackoverflow.com/a/23212515/717419
    if platform.system() == 'Windows':
        raise SkipTest
    indexed = []
    with open('data/genes.fasta') as genes:
        fasta = genes.readlines()
    n_lines = sum(1 for line in fasta)
    for n in range(n_lines):
        with NamedTemporaryFile(mode='w') as lines:
            for i, line in enumerate(fasta):
                if i == n and line[0] != '>' and len(line) == 71:
                    line = line[:-3] + '\n'
                    full_line = True
                elif i == n:
                    full_line = False
                lines.write(line)
                lines.flush()
            name = lines.name
            if full_line:
                try:
                    Faidx(name)
                    indexed.append(True)
                except FastaIndexingError:
                    indexed.append(False)
    assert not any(indexed)

def test_valgrind_long_lines(remove_index):
    """ Makes all full-length lines long and checks that error is raised
    in all appropriate circumstances.
    """
    # http://stackoverflow.com/a/23212515/717419
    if platform.system() == 'Windows':
        raise SkipTest
    indexed = []
    with open('data/genes.fasta') as genes:
        fasta = genes.readlines()
    n_lines = sum(1 for line in fasta)
    for n in range(n_lines):
        with NamedTemporaryFile(mode='w') as lines:
            for i, line in enumerate(fasta):
                if i == n and line[0] != '>' and len(line) == 71:
                    line = line.rstrip('\n') + 'NNN' + '\n'
                    full_line = True
                elif i == n:
                    full_line = False
                lines.write(line)
                lines.flush()
            name = lines.name
            if full_line:
                try:
                    Faidx(name)
                    indexed.append(True)
                except FastaIndexingError:
                    indexed.append(False)
    assert not any(indexed)

def test_valgrind_blank_lines(remove_index):
    """ Makes all full-length lines blank and checks that error is raised
    in all appropriate circumstances.
    """
    # http://stackoverflow.com/a/23212515/717419
    if platform.system() == 'Windows':
        raise SkipTest
    indexed = []
    with open('data/genes.fasta') as genes:
        fasta = genes.readlines()
    n_lines = sum(1 for line in fasta)
    for n in range(n_lines):
        with NamedTemporaryFile(mode='w') as lines:
            for i, line in enumerate(fasta):
                if i == n and line[0] != '>' and len(line) == 71:
                    line = '\n'
                    full_line = True
                elif i == n:
                    full_line = False
                lines.write(line)
                lines.flush()
            name = lines.name
            if full_line:
                try:
                    Faidx(name)
                    indexed.append(True)
                except FastaIndexingError:
                    indexed.append(False)
    assert not any(indexed)

def test_reindex_on_modification(remove_index):
    """ This test ensures that the index is regenerated when the FASTA
    modification time is newer than the index modification time.
    mdshw5/pyfaidx#50 """
    faidx = Faidx('data/genes.fasta')
    index_mtime = getmtime(faidx.indexname)
    faidx.close()
    os.utime('data/genes.fasta', (index_mtime + 10, ) * 2)
    time.sleep(2)
    faidx = Faidx('data/genes.fasta')
    assert getmtime(faidx.indexname) > index_mtime

def test_build_issue_83(remove_index):
    """ Ensure that blank lines between entries are treated in the
    same way as samtools 1.2. See mdshw5/pyfaidx#83.
    """
    expect_index = ("MT	119	4	70	71\nGL000207.1	60	187	60	61\n")
    index_file = Faidx('data/issue_83.fasta').indexname
    result_index = open(index_file).read()
    os.remove('data/issue_83.fasta.fai')
    assert result_index == expect_index

def test_build_issue_96_fail_build_faidx(remove_index):
    """ Ensure that the fasta file is closed if construction of the 'Faidx' file
    when attempting to build an index.
    See mdshw5/pyfaidx#96
    """
    tmp_dir = mkdtemp()
    try:
        fasta_path = os.path.join(tmp_dir, 'issue_96.fasta')
        # Write simple fasta file with inconsistent sequence line lengths,
        # so building an index raises a 'FastaIndexingError'
        with open(fasta_path, 'w') as fasta_out:
            fasta_out.write(
                ">seq1\nCTCCGGGCCCAT\nAACACTTGGGGGTAGCTAAAGTGAA\nATAAAGCCTAAA\n"
            )

        builtins_open = builtins.open

        opened_files = []

        def test_open(*args, **kwargs):
            f = builtins_open(*args, **kwargs)
            opened_files.append(f)
            return f

        with mock.patch('builtins.open', side_effect=test_open):
            try:
                Faidx(fasta_path)
                remove_index.assertFail(
                    "Faidx construction should fail with 'FastaIndexingError'."
                )
            except FastaIndexingError:
                pass
        assert all(f.closed for f in opened_files)
    finally:
        shutil.rmtree(tmp_dir)

def test_build_issue_96_fail_read_malformed_index_duplicate_key(remove_index):
    """ Ensure that the fasta file is closed if construction of the 'Faidx' file
    fails when attempting to read a pre-existing index. The index is malformed because
    it contains mulitple occurrences of the same index.
    See mdshw5/pyfaidx#96
    """
    tmp_dir = mkdtemp()
    try:
        fasta_path = os.path.join(tmp_dir, 'issue_96.fasta')
        faidx_path = os.path.join(tmp_dir, 'issue_96.fasta.fai')
        # Write simple fasta file
        with open(fasta_path, 'w') as fasta_out:
            fasta_out.write(">seq1\nCTCCGGGCCCAT\nATAAAGCCTAAA\n")
        with open(faidx_path, 'w') as faidx_out:
            faidx_out.write("seq1\t24\t6\t12\t13\nseq1\t24\t6\t12\t13\n")

        builtins_open = builtins.open

        opened_files = []

        def test_open(*args, **kwargs):
            f = builtins_open(*args, **kwargs)
            opened_files.append(f)
            return f

        with mock.patch('builtins.open', side_effect=test_open):
            try:
                Faidx(fasta_path)
                remove_index.assertFail(
                    "Faidx construction should fail with 'ValueError'.")
            except ValueError:
                pass
        assert all(f.closed for f in opened_files)
    finally:
        shutil.rmtree(tmp_dir)

def test_read_back_index(remove_index):
    """Ensure that index files written with write_fai() can be read back"""
    import locale
    import platform
    
    if platform.system() == "Linux":
        new_locale = 'en_US.utf8'
    elif platform.system() == "Darwin":
        new_locale = 'en_US.UTF-8'
    
    old_locale = locale.getlocale(locale.LC_NUMERIC)
    try:
        locale.setlocale(locale.LC_NUMERIC, new_locale)
        faidx = Faidx('data/genes.fasta')
        faidx.write_fai()
        faidx = Faidx('data/genes.fasta', build_index=False)
    finally:
        locale.setlocale(locale.LC_NUMERIC, old_locale)

def test_issue_134_no_build_index(remove_index):
    """ Ensure that index file is not built when build_index=False. See mdshw5/pyfaidx#134.
    """
    with pytest.raises(IndexNotFoundError):
        faidx = Faidx('data/genes.fasta', build_index=False)

def test_issue_144_no_defline(remove_index):
    """ Ensure that an exception is raised when a file contains no deflines. See mdshw5/pyfaidx#144.
    """
    tmp_dir = mkdtemp()
    try:
        fasta_path = os.path.join(tmp_dir, 'issue_144.fasta')
        # Write simple fasta file
        with open(fasta_path, 'w') as fasta_out:
            fasta_out.write("CTCCGGGCCCAT\nATAAAGCCTAAA\n")
        with pytest.raises(FastaIndexingError):
            faidx = Faidx(fasta_path)
    finally:
        shutil.rmtree(tmp_dir)

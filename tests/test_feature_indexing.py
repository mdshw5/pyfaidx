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
    # Removed old GI-based expect_index assignments
    expect_index = (
        "AB821309.1\t3510\t96\t70\t71\n"
        "KF435150.1\t481\t3754\t70\t71\n"
        "KF435149.1\t642\t4316\t70\t71\n"
        "NR_104216.1\t4573\t5071\t70\t71\n"
        "NR_104215.1\t5317\t9813\t70\t71\n"
        "NR_104214.1\t1869\t15327\t70\t71\n"
        "NR_104208.1\t1414\t17330\t70\t71\n"
        "NR_104206.1\t1414\t18853\t70\t71\n"
        "NR_104205.1\t1298\t20376\t70\t71\n"
        "NR_104160.1\t923\t21795\t70\t71\n"
        "NR_104158.1\t399\t22839\t70\t71\n"
        "XR_241081.1\t1009\t23362\t70\t71\n"
        "XR_241080.1\t4884\t24504\t70\t71\n"
        "XR_241079.1\t2819\t29576\t70\t71\n"
        "XR_241078.1\t567\t32553\t70\t71\n"
        "XR_241065.1\t2806\t33226\t70\t71\n"
        "XR_241064.1\t2824\t36170\t70\t71\n"
        "XR_241055.1\t4111\t39123\t70\t71\n"
        "XR_241054.1\t3251\t43384\t70\t71\n"
        "XR_241053.1\t909\t46770\t70\t71\n"
    )
    index_file = Faidx('data/genes.fasta').indexname
    result_index = open(index_file).read()
    assert result_index == expect_index

def test_build_issue_141(remove_index):
    expect_index = (
        "AB821309.1\t3510\t97\t70\t72\n"
        "KF435150.1\t481\t3807\t70\t72\n"
        "KF435149.1\t642\t4377\t70\t72\n"
        "NR_104216.1\t4573\t5143\t70\t72\n"
        "NR_104215.1\t5317\t9952\t70\t72\n"
        "NR_104214.1\t1869\t15543\t70\t72\n"
        "NR_104208.1\t1414\t17574\t70\t72\n"
        "NR_104206.1\t1414\t19119\t70\t72\n"
        "NR_104205.1\t1298\t20664\t70\t72\n"
        "NR_104160.1\t923\t22103\t70\t72\n"
        "NR_104158.1\t399\t23162\t70\t72\n"
        "XR_241081.1\t1009\t23692\t70\t72\n"
        "XR_241080.1\t4884\t24850\t70\t72\n"
        "XR_241079.1\t2819\t29993\t70\t72\n"
        "XR_241078.1\t567\t33012\t70\t72\n"
        "XR_241065.1\t2806\t33695\t70\t72\n"
        "XR_241064.1\t2824\t36681\t70\t72\n"
        "XR_241055.1\t4111\t39676\t70\t72\n"
        "XR_241054.1\t3251\t43997\t70\t72\n"
        "XR_241053.1\t909\t47431\t70\t72\n"
    )
    index_file = Faidx('data/issue_141.fasta').indexname
    result_index = open(index_file).read()
    os.remove('data/issue_141.fasta.fai')
    print(result_index)
    assert result_index == expect_index

def test_build_issue_111(remove_index):
    # Removed old GI-based expect_index assignments
    expect_index = (
        "AB821309\t3510\t96\t70\t71\n"
        "KF435150\t481\t3754\t70\t71\n"
        "KF435149\t642\t4316\t70\t71\n"
        "NR_104216\t4573\t5071\t70\t71\n"
        "NR_104215\t5317\t9813\t70\t71\n"
        "NR_104214\t1869\t15327\t70\t71\n"
        "NR_104208\t1414\t17330\t70\t71\n"
        "NR_104206\t1414\t18853\t70\t71\n"
        "NR_104205\t1298\t20376\t70\t71\n"
        "NR_104160\t923\t21795\t70\t71\n"
        "NR_104158\t399\t22839\t70\t71\n"
        "XR_241081\t1009\t23362\t70\t71\n"
        "XR_241080\t4884\t24504\t70\t71\n"
        "XR_241079\t2819\t29576\t70\t71\n"
        "XR_241078\t567\t32553\t70\t71\n"
        "XR_241065\t2806\t33226\t70\t71\n"
        "XR_241064\t2824\t36170\t70\t71\n"
        "XR_241055\t4111\t39123\t70\t71\n"
        "XR_241054\t3251\t43384\t70\t71\n"
        "XR_241053\t909\t46770\t70\t71\n"
    )
    index = Faidx(
        'data/genes.fasta',
        read_long_names=True,
        key_function=lambda x: x.split('.')[0])
    result_index = ''.join(index._index_as_string())
    assert result_index == expect_index

def test_order(remove_index):
    # Removed old GI-based order assignment
    order = (
       'AB821309.1', 'KF435150.1', 'KF435149.1', 'NR_104216.1', 'NR_104215.1', 'NR_104214.1', 'NR_104208.1', 'NR_104206.1', 'NR_104205.1', 'NR_104160.1', 'NR_104158.1', 'XR_241081.1', 'XR_241080.1', 'XR_241079.1', 'XR_241078.1', 'XR_241065.1', 'XR_241064.1', 'XR_241055.1', 'XR_241054.1', 'XR_241053.1'
    )
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

@pytest.mark.xfail(raises=IndexNotFoundError)
def test_issue_134_no_build_index(remove_index):
    """ Ensure that index file is not built when build_index=False. See mdshw5/pyfaidx#134.
    """
    faidx = Faidx('data/genes.fasta', build_index=False)

@pytest.mark.xfail(raises=FastaIndexingError)
def test_issue_144_no_defline():
    """ Ensure that an exception is raised when a file contains no deflines. See mdshw5/pyfaidx#144.
    """
    tmp_dir = mkdtemp()
    try:
        fasta_path = os.path.join(tmp_dir, 'issue_144.fasta')
        # Write simple fasta file
        with open(fasta_path, 'w') as fasta_out:
            fasta_out.write("CTCCGGGCCCAT\nATAAAGCCTAAA\n")
        faidx = Faidx(fasta_path)
    finally:
        shutil.rmtree(tmp_dir)

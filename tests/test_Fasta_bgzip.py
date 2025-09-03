import os
import pytest
from pyfaidx import Fasta, Faidx, UnsupportedCompressionFormat, FetchError
from itertools import chain

path = os.path.dirname(__file__)
os.chdir(path)

try:
    from Bio import SeqIO
    bio = True
except ImportError:
    bio = False
    
@pytest.fixture
def remove_index():
    yield
    try:
        os.remove('data/genes.fasta.gz.fai')
        os.remove('data/genes.fasta.gz.gzi')
    except EnvironmentError:
        pass  # some tests may delete this file

@pytest.mark.skipif(not bio, reason="Biopython is not installed.")
def test_build_issue_126():
    """ Samtools BGZF index should be identical to pyfaidx BGZF index """
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
    index_file = Faidx('data/genes.fasta.gz').indexname
    result_index = open(index_file).read()
    assert result_index == expect_index

def test_integer_slice():
    fasta = Fasta('data/genes.fasta.gz')
    expect = fasta['AB821309.1'][:100].seq
    result = fasta[0][:100].seq
    assert expect == result

def test_integer_index():
    fasta = Fasta('data/genes.fasta.gz')
    expect = fasta['AB821309.1'][100].seq
    result = fasta[0][100].seq
    assert expect == result

def test_fetch_whole_fasta():
    expect = [line.rstrip('\n') for line in open('data/genes.fasta') if line[0] != '>']
    result = list(chain(*([line for line in record] for record in Fasta('data/genes.fasta.gz', as_raw=True))))
    assert expect == result

def test_line_len():
    fasta = Fasta('data/genes.fasta.gz')
    for record in fasta:
        assert len(next(iter(record))) == fasta.faidx.index[record.name].lenc

@pytest.mark.xfail(raises=UnsupportedCompressionFormat)
def test_mutable_bgzf():
    fasta = Fasta('data/genes.fasta.gz', mutable=True)

def test_long_names():
    """ Test that deflines extracted using FastaRecord.long_name are
    identical to deflines in the actual file.
    """
    deflines = []
    with open('data/genes.fasta') as fasta_file:
        for line in fasta_file:
            if line[0] == '>':
                deflines.append(line[1:-1])
    fasta = Fasta('data/genes.fasta.gz')
    long_names = []
    for record in fasta:
        long_names.append(record.long_name)
    assert deflines == long_names

def test_fetch_whole_entry():
    faidx = Faidx('data/genes.fasta.gz')
    expect = ('ATGACATCATTTTCCACCTCTGCTCAGTGTTCAACATCTGA'
            'CAGTGCTTGCAGGATCTCTCCTGGACAAATCAATCAGGTACGACCA'
            'AAACTGCCGCTTTTGAAGATTTTGCATGCAGCAGGTGCGCAAGG'
            'TGAAATGTTCACTGTTAAAGAGGTCATGCACTATTTAGGTCAGTACAT'
            'AATGGTGAAGCAACTTTATGATCAGCAGGAGCAGCATATGGTATATTG'
            'TGGTGGAGATCTTTTGGGAGAACTACTGGGACGTCAGAGCTTCTCCGTG'
            'AAAGACCCAAGCCCTCTCTATGATATGCTAAGAAAGAATCTTGTCACTTT'
            'AGCCACTGCTACTACAGCAAAGTGCAGAGGAAAGTTCCACTTCCAGAAAAA'
            'GAACTACAGAAGACGATATCCCCACACTGCCTACCTCAGAGCATAAATGCA'
            'TACATTCTAGAGAAGGTGATTGAAGTGGGAAAAAATGATGACCTGGAGGACTC')
    result = faidx.fetch('KF435150.1',
                         1, 481)
    assert str(result) == expect

def test_fetch_middle():
    faidx = Faidx('data/genes.fasta.gz')
    expect = 'TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA'
    result = faidx.fetch('KF435150.1',
                         100, 150)
    assert str(result) == expect

def test_fetch_end():
    faidx = Faidx('data/genes.fasta.gz')
    expect = 'TC'
    result = faidx.fetch('KF435150.1',
                         480, 481)
    assert str(result) == expect

@pytest.mark.xfail(raises=FetchError)
def test_fetch_border():
    """ Fetch past the end of a gene entry """
    faidx = Faidx('data/genes.fasta.gz')
    expect = 'TC'
    result = faidx.fetch('KF435150.1', 480, 500)

def test_rev():
    faidx = Faidx('data/genes.fasta.gz')
    expect = 'GA'
    result = faidx.fetch('KF435150.1',
                         480, 481)
    assert str(-result) == expect, result

@pytest.mark.xfail(raises=FetchError)
def test_fetch_past_bounds():
    """ Fetch past the end of a gene entry """
    faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
    result = faidx.fetch('KF435150.1', 480, 5000)

@pytest.mark.xfail(raises=FetchError)
def test_fetch_negative():
    """ Fetch starting with a negative coordinate """
    faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
    result = faidx.fetch('KF435150.1', -10, 10)

@pytest.mark.xfail(raises=FetchError)
def test_fetch_reversed_coordinates():
    """ Fetch starting with a negative coordinate """
    faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
    result = faidx.fetch('KF435150.1', 50, 10)

@pytest.mark.xfail(raises=FetchError)
def test_fetch_keyerror():
    """ Fetch a key that does not exist """
    faidx = Faidx('data/genes.fasta.gz', strict_bounds=True)
    result = faidx.fetch('KFXXXXXXX.1|', 1, 10)

def test_blank_string():
    """ seq[0:0] should return a blank string mdshw5/pyfaidx#53 """
    fasta = Fasta('data/genes.fasta.gz', as_raw=True)
    assert fasta['KF435150.1'][0:0] == ''

def test_slice_from_beginning():
    fasta = Fasta('data/genes.fasta.gz', as_raw=True)
    assert fasta['KF435150.1'][:4] == 'ATGA'

def test_slice_from_end():
    fasta = Fasta('data/genes.fasta.gz', as_raw=True)
    assert fasta['KF435150.1'][-4:] == 'ACTC'

def test_issue_74_start():
    f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
    f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
    assert f0['KF435150.1'][0:90].start == f1['KF435150.1'][0:90].start - 1

def test_issue_74_consistency():
    f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
    f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
    assert str(f0['KF435150.1'][0:90]) == str(f1['KF435150.1'][0:90])

def test_issue_74_end_faidx():
    f0 = Faidx('data/genes.fasta.gz', one_based_attributes=False)
    f1 = Faidx('data/genes.fasta.gz', one_based_attributes=True)
    end0 = f0.fetch('KF435150.1', 1, 90).end
    end1 = f1.fetch('KF435150.1', 1, 90).end
    assert end0 == end1

def test_issue_74_end_fasta():
    f0 = Fasta('data/genes.fasta.gz', one_based_attributes=False)
    f1 = Fasta('data/genes.fasta.gz', one_based_attributes=True)
    end0 = f0['KF435150.1'][1:90].end
    end1 = f1['KF435150.1'][1:90].end
    print((end0, end1))
    assert end0 == end1

def test_issue_79_fix():
    f = Fasta('data/genes.fasta.gz')
    s = f['KF435150.1'][100:105]
    print((s.start, s.end))
    assert (101, 105) == (s.start, s.end)

def test_issue_79_fix_negate():
    f = Fasta('data/genes.fasta.gz')
    s = f['KF435150.1'][100:105]
    s = -s
    print((s.start, s.end))
    assert (105, 101) == (s.start, s.end)

def test_issue_79_fix_one_based_false():
    f = Fasta('data/genes.fasta.gz', one_based_attributes=False)
    s = f['KF435150.1'][100:105]
    print((s.start, s.end))
    assert (100, 105) == (s.start, s.end)

def test_issue_79_fix_one_based_false_negate():
    f = Fasta('data/genes.fasta.gz', one_based_attributes=False)
    s = f['KF435150.1'][100:105]
    print(s.__dict__)
    s = -s
    print(s.__dict__)
    assert (105, 100) == (s.start, s.end)

@pytest.mark.xfail(raises=FetchError)
def test_fetch_border_padded(remove_index):
    """ Fetch past the end of a gene entry """
    faidx = Faidx('data/genes.fasta.gz', default_seq='N')
    expect = 'TCNNNNNNNNNNNNNNNNNNN'
    result = faidx.fetch('KF435150.1',
                             480, 500)
    
@pytest.mark.xfail(raises=UnsupportedCompressionFormat)
def test_mutable_bgzf_fasta(remove_index):
    """ Test that mutable Fasta raises error with bgzf files """
    fasta = Fasta('data/genes.fasta.gz', mutable=True)

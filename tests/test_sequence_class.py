import pytest
from pyfaidx import Sequence, complement

seq = Sequence(name='KF435150.1', seq='TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA',
                    start=100, end=150)

seq_invalid = Sequence(name='KF435150.1', seq='TTGAAGATTTPGCATGCAGCAGGTGCGCAAGGTGAAATNTTCACTGTTAAA',
                    start=100, end=150)

comp_valid = 'TTGAAGATTTnGCATGCAGCAGGtgccaAGGTGAAATGTTNACTGTTAAA'

comp_invalid = 'TTGAAGATTTnGCATGCAGCPQGtgccaAGGTGAAATGTTNACTGTTAAA'

def test_negate():
    assert str(-seq) == str(seq.complement[::-1])

def test_negate_metadata():
    # Negate should affect __repr__ the same way as reverse and complement
    seq_neg = -seq
    assert seq_neg.__repr__() == seq.complement[::-1].__repr__()

@pytest.mark.xfail(raises=ValueError)
def test_seq_invalid():
    seq_invalid.complement()

def test_integer_index():
    assert seq[1].seq == 'T'

def test_slice_index():
    assert seq[0:10].seq == 'TTGAAGATTT'

@pytest.mark.xfail(raises=ValueError)
def test_comp_invalid():
    complement(comp_invalid)

@pytest.mark.xfail(raises=ValueError)
def test_check_coordinates():
    x = Sequence(name='KF435150.1', seq='TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA',
                 start=100, end=110)
    _ = x[:]

def test_comp_valid():
    assert complement(comp_valid).startswith("AACTTCTAAAnCG")
    assert complement(complement(comp_valid)) == comp_valid

def test_comp_empty():
    assert complement('') == ''

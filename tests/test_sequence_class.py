from pyfaidx import Sequence, complement
from nose.tools import assert_raises, raises

seq = Sequence(name='gi|557361099|gb|KF435150.1|', seq='TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA',
                    start=100, end=150)

seq_invalid = Sequence(name='gi|557361099|gb|KF435150.1|', seq='TTGAAGATTTPGCATGCAGCAGGTGCGCAAGGTGAAATNTTCACTGTTAAA',
                    start=100, end=150)
                    
seq_0 = Sequence(name='gi|557361099|gb|KF435150.1|', seq='TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA',
                    start=100, end=150, one_based_attributes=False)                    

comp_valid = 'TTGAAGATTTnGCATGCAGCAGGtgccaAGGTGAAATGTTNACTGTTAAA'

comp_invalid = 'TTGAAGATTTnGCATGCAGCPQGtgccaAGGTGAAATGTTNACTGTTAAA'

def test_negate():
    assert str(-seq) == str(seq.complement[::-1])

def test_negate_metadata():
    # Negate should affect __repr__ the same way as reverse and complement
    seq_neg = -seq
    assert seq_neg.__repr__() == seq.complement[::-1].__repr__()

def test_seq_invalid():
    assert_raises(ValueError, lambda: seq_invalid.complement)

def test_integer_index():
    assert seq[1].seq == 'T'

def test_slice_index():
    assert seq[0:10].seq == 'TTGAAGATTT'

@raises(ValueError)
def test_comp_invalid():
    complement(comp_invalid)

def test_comp_valid():
    complement(comp_valid)
    
def test_zero_start():
    assert seq_0.start == seq.start - 1

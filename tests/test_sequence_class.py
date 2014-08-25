from pyfaidx import Sequence

test_seq = Sequence(name='KF435150.1', seq='TTGAAGATTTTGCATGCAGCAGGTGCGCAAGGTGAAATGTTCACTGTTAAA',
                    start=100, end=150)

def test_negate():
    assert str(-test_seq) == str(test_seq.complement[::-1])

def test_negate_metadata():
    # Negate should affect __repr__ the same way as reverse and complement
    test_seq_neg = -test_seq
    assert test_seq_neg.__repr__() == test_seq.complement[::-1].__repr__()

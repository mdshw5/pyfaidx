import os
import tempfile
import pytest
from pyfaidx import wrap_sequence

def test_wrap_sequence_preserves_newline():
    seq = 'ACTG' * 5  # 20 bases
    # Try all common newlines
    for newline in ['\n', '\r\n', '\r']:
        lines = list(wrap_sequence(7, seq, newline=newline))
        # All lines except last should end with the correct newline
        for l in lines:
            assert l.endswith(newline)
        # Join and split to check round-trip
        joined = ''.join(lines)
        split = joined.split(newline)
        # The last split element should be empty (trailing newline)
        assert split[-1] == ''
        # The sequence reconstructed (without newlines) should match original
        assert ''.join(split[:-1]) == seq

if __name__ == '__main__':
    pytest.main([__file__])

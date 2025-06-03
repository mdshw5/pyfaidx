import os
import tempfile
import pytest
from pyfaidx import Fasta

def make_long_fasta(tmp_path, seq_len=6_000_000):
    fasta_path = tmp_path / "longseq.fa"
    with open(fasta_path, "w") as f:
        f.write(">chr1\n")
        # Write sequence in 60bp lines (standard wrapping)
        seq = ("ACGT" * (seq_len // 4 + 1))[:seq_len]
        for i in range(0, seq_len, 60):
            f.write(seq[i:i+60] + "\n")
    return fasta_path

@pytest.mark.parametrize("repeat", range(5))
def test_long_fasta_entry_length(tmp_path, repeat):
    fasta_path = make_long_fasta(tmp_path)
    fa = Fasta(str(fasta_path))
    # Read by name
    seq1 = fa["chr1"][:].seq.upper()
    assert len(seq1) == 6_000_000
    # Read by explicit range
    seq2 = fa["chr1"][0:6_000_000].seq.upper()
    assert len(seq2) == 6_000_000
    # Read by index (should be only one entry)
    seq3 = fa[0][:].seq.upper()
    assert len(seq3) == 6_000_000
    # All reads should be identical
    assert seq1 == seq2 == seq3

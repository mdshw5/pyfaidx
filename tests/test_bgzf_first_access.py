import pytest
from pathlib import Path
from pyfaidx import Fasta

TEST_FASTA = Path(__file__).parent / "data/Anomala_cuprea_entomopoxvirus.faa.gz"

@pytest.fixture(autouse=True)
def teardown():
    # Setup: ensure the test fasta and its index are present
    yield
    # Teardown: remove the index files after each test
    for ext in ['.fai', '.gzi']:
        index_file = TEST_FASTA.with_suffix(TEST_FASTA.suffix + ext)
        if index_file.exists():
            index_file.unlink()

@pytest.mark.parametrize("seq_id, slc, expected", [
    ("YP_009001474.1", slice(10, 20), "DKKDKELYIM"),
])
def test_bgzf_first_access_valueerror(seq_id, slc, expected):
    # First access: should not raise ValueError
    fasta = Fasta(TEST_FASTA)
    try:
        seq = fasta[seq_id][slc.start:slc.stop]
    except ValueError as e:
        pytest.fail(f"Unexpected ValueError on first access: {e}")
    assert str(seq) == expected

    # Second access: should also work
    fasta = Fasta(TEST_FASTA)
    seq = fasta[seq_id][slc.start:slc.stop]
    assert str(seq) == expected

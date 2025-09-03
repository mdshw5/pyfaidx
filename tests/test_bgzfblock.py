import struct
import pytest
from pyfaidx import BgzfBlock, make_virtual_offset

def test_bgzfblock_fields():
    block = BgzfBlock(10, 20, 30, 40)
    assert block.cstart == 10
    assert block.clen == 20
    assert block.ustart == 30
    assert block.ulen == 40

def test_bgzfblock_getitem():
    block = BgzfBlock(1, 2, 3, 4)
    assert block[0] == 1
    assert block[1] == 2
    assert block[2] == 3
    assert block[3] == 4
    assert block['cstart'] == 1
    assert block['clen'] == 2
    assert block['ustart'] == 3
    assert block['ulen'] == 4

def test_bgzfblock_as_bytes():
    block = BgzfBlock(123, 456, 789, 1011)
    b = block.as_bytes()
    # Only cstart and ustart are packed
    expected = struct.pack('<QQ', 123, 789)
    assert b == expected

def test_bgzfblock_lt():
    block = BgzfBlock(0, 0, 100, 10)
    assert (block < 200) is True
    assert (block < 50) is False

def test_bgzfblock_len():
    block = BgzfBlock(0, 0, 0, 42)
    assert len(block) == 42

def test_bgzfblock_empty_true():
    block = BgzfBlock(0, 0, 0, 0)
    assert block.empty is True

def test_bgzfblock_empty_false():
    block = BgzfBlock(0, 0, 0, 1)
    assert block.empty is False

def test_make_virtual_offset():
    # Normal case
    assert make_virtual_offset(1, 1) == (1 << 16) | 1
    assert make_virtual_offset(0, 0) == 0
    assert make_virtual_offset(12345, 54321) == (12345 << 16) | 54321
    # Edge cases
    assert make_virtual_offset(0, 65535) == (0 << 16) | 65535
    assert make_virtual_offset(281474976710655, 0) == (281474976710655 << 16)
    # Error cases
    import pytest
    with pytest.raises(ValueError):
        make_virtual_offset(-1, 0)
    with pytest.raises(ValueError):
        make_virtual_offset(0, -1)
    with pytest.raises(ValueError):
        make_virtual_offset(0, 65536)
    with pytest.raises(ValueError):
        make_virtual_offset(281474976710656, 0)

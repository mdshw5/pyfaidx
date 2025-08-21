import struct
import pytest
from pyfaidx import BgzfBlock

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

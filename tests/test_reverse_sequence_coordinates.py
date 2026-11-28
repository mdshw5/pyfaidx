import pytest

from pyfaidx import Sequence


@pytest.mark.parametrize('one_based', [False, True])
@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('comp', [False, True])
@pytest.mark.parametrize('region', [
    slice(8, 2, -1), slice(5, 1, -1), slice(7, None, -1),
    slice(None, 3, -1), slice(-2, -7, -1), slice(None, None, -1),
])
def test_reverse_slice_coordinates(one_based, reverse, comp, region):
    text = 'ACGTTAGCCTAA'
    start = 100
    end = start + len(text) - int(one_based)
    sequence = Sequence('record', text, start, end)
    # Track each base's endpoint labels independently of Sequence arithmetic.
    labels = [(i, i + int(not one_based))
              for i in range(start, start + len(text))]
    if reverse:
        sequence = sequence.reverse
        labels = [pair[::-1] for pair in labels[::-1]]
    if comp:
        sequence = sequence.complement

    result = sequence[region]
    selected = [pair[::-1] for pair in labels[region]]
    assert result.seq == sequence.seq[region]
    assert (result.start, result.end) == (selected[0][0], selected[-1][1])
    assert result.comp == comp
    assert result[:].seq == result.seq
    assert (result[:].start, result[:].end) == (result.start, result.end)
    assert result.reverse.reverse.fancy_name == result.fancy_name

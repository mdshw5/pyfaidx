import pytest

from pyfaidx import Sequence


@pytest.mark.parametrize('one_based', [False, True])
@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('comp', [False, True])
@pytest.mark.parametrize('index', [0, 2, -1, -5])
def test_integer_coordinates_match_selected_base(one_based, reverse, comp, index):
    start = 1 if one_based else 0
    sequence = Sequence('chr1', 'ACGTA', start=start, end=5)
    labels = [(i, i + int(not one_based)) for i in range(start, start + 5)]
    if reverse:
        sequence = sequence.reverse
        labels = [pair[::-1] for pair in labels[::-1]]
    if comp:
        sequence = sequence.complement
    selected = sequence[index]
    assert selected.seq == sequence.seq[index]
    assert (selected.start, selected.end) == labels[index]
    assert selected.comp is comp
    assert selected[:].fancy_name == selected.fancy_name


@pytest.mark.parametrize('comp', [False, True])
def test_integer_without_coordinates_preserves_complement(comp):
    sequence = Sequence('chr1', 'ACGTA', comp=comp)
    selected = sequence[-1]
    assert selected.seq == 'A'
    assert selected.start is None
    assert selected.end is None
    assert selected.comp is comp


@pytest.mark.parametrize('index', [-6, -10, 5, 10])
def test_integer_outside_sequence_raises(index):
    sequence = Sequence('chr1', 'ACGTA', start=1, end=5)
    with pytest.raises(IndexError):
        sequence[index]


@pytest.mark.parametrize('index', [-1, 0])
def test_integer_empty_sequence_raises(index):
    with pytest.raises(IndexError):
        Sequence('chr1', '')[index]

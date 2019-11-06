from pfsspy import FieldLine
import pytest


@pytest.mark.parametrize('x, open, pol',
                         [[[1, 2.5], True, 1],
                          [[2.5, 1], True, -1],
                          [[1, 1], False, 0],
                          ])
def test_open(x, open, pol):
    fline = FieldLine(x, [0, 0], [0, 0], None, None)

    assert (fline.is_open == open)
    assert (fline.polarity == pol)

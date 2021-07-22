from pfsspy.fieldline import FieldLine, FieldLines, OpenFieldLines, ClosedFieldLines


import pytest


@pytest.mark.parametrize('x, open, pol',
                         [[[1, 2.5], True, 1],
                          [[2.5, 1], True, -1],
                          [[1, 1], False, 0],
                          ])
def test_open(x, open, pol):
    fline = FieldLine(x, [0, 0], [0, 0], None)

    assert (fline.is_open == open)
    assert (fline.polarity == pol)
    assert len(fline) == 2

    flines = FieldLines([fline])

    assert len(flines.open_field_lines) == int(open)
    assert len(flines.closed_field_lines) == int(not open)


@pytest.mark.parametrize('x, cls',
                         [[[1, 2.5], ClosedFieldLines],
                          [[1, 1], OpenFieldLines],
                          ])
def test_flines_errors(x, cls):
    fline = FieldLine(x, [0, 0], [0, 0], None)
    with pytest.raises(ValueError):
        cls([fline])

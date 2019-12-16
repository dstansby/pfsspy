from pfsspy.fieldline import FieldLine, FieldLines, OpenFieldLines, ClosedFieldLines

import astropy.units as u
from astropy.time import Time
from sunpy.coordinates import frames as sunframes

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

    flines = FieldLines([fline])

    assert len(flines.open_field_lines) == int(open)
    assert len(flines.closed_field_lines) == int(not open)


@pytest.mark.parametrize('x, cls',
                         [[[1, 2.5], ClosedFieldLines],
                          [[1, 1], OpenFieldLines],
                          ])
def test_flines_errors(x, cls):
    fline = FieldLine(x, [0, 0], [0, 0], None, None)
    with pytest.raises(ValueError):
        cls([fline])


def test_transform():
    # Check that field lines can be transformed into different coordinates
    obstime = Time('1992-12-21')
    stonyhurst = sunframes.HeliographicStonyhurst(
        12 * u.deg, 0 * u.deg, obstime=obstime)
    fline = FieldLine([1, 2.5], [0, 0], [0, 0], None, None)
    # Check field line transform
    fline.coords.transform_to(stonyhurst)

    # Check footpoint transform
    fline.solar_footpoint.transform_to(stonyhurst)

from pfsspy import FieldLine, FieldLines

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


def test_fieldlines_obstime():
    # Check that obstime is propagated to FieldLines and subsequent preoperties
    x = [1, 2.5]
    obstime = Time('1992-12-21')
    fline1 = FieldLine(x, [0, 0], [0, 0], obstime, None)
    fline2 = FieldLine(x, [0, 0], [0, 0], obstime, None)

    flines = FieldLines([fline1, fline2])

    assert flines.solar_feet.obstime == obstime
    assert flines.source_surface_feet.obstime == obstime


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

from pfsspy.fieldline import FieldLine, FieldLines, OpenFieldLines, ClosedFieldLines

import astropy.units as u
import astropy.constants as const
from astropy.time import Time
from astropy.coordinates import SkyCoord
from sunpy.coordinates import frames as sunframes

import pytest

rsun = const.R_sun


@pytest.mark.parametrize('r, open, pol',
                         [[[1, 2.5], True, 1],
                          [[2.5, 1], True, -1],
                          [[1, 1], False, 0],
                          ])
def test_open(r, open, pol):
    coord = SkyCoord(0 * u.deg, 0 * u.deg, r * rsun,
                     frame='heliographic_carrington')
    fline = FieldLine(coord, None)

    assert (fline.is_open == open)
    assert (fline.polarity == pol)

    flines = FieldLines([fline])

    assert len(flines.open_field_lines) == int(open)
    assert len(flines.closed_field_lines) == int(not open)


@pytest.mark.parametrize('r, cls',
                         [[[1, 2.5], ClosedFieldLines],
                          [[1, 1], OpenFieldLines],
                          ])
def test_flines_errors(r, cls):
    coord = SkyCoord(0 * u.deg, 0 * u.deg, r * rsun,
                     frame='heliographic_carrington')
    fline = FieldLine(coord, None)
    with pytest.raises(ValueError):
        cls([fline])


def test_transform():
    # Check that field lines can be transformed into different coordinates
    obstime = Time('1992-12-21')
    stonyhurst = sunframes.HeliographicStonyhurst(
        12 * u.deg, 0 * u.deg, obstime=obstime)
    coord = SkyCoord([0 * u.deg, 0 * u.deg],
                     [0 * u.deg, 0 * u.deg],
                     [1 * rsun, 2.5 * rsun],
                     frame='heliographic_carrington',
                     observer='earth',
                     obstime=obstime)
    fline = FieldLine(coord, None)
    # Check field line transform
    fline.coords.transform_to(stonyhurst)

    # Check footpoint transform
    fline.solar_footpoint.transform_to(stonyhurst)

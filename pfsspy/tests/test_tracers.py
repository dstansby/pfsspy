import astropy.coordinates as coord
import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest

import pfsspy
from pfsspy import tracing

from .example_maps import dipole_map, dipole_result


@pytest.fixture(params=[tracing.PythonTracer(),
                        tracing.FortranTracer()])
def flines(dipole_result, request):
    tracer = request.param
    input, out = dipole_result
    out_frame = out.coordinate_frame

    seed = coord.SkyCoord(2*u.deg, -45*u.deg, 1.01*const.R_sun, frame=out_frame)
    flines = tracer.trace(seed, out)
    print(flines[0].coords)
    return flines


def test_field_lines(flines):
    assert isinstance(flines[0],
                      pfsspy.fieldline.FieldLine)
    assert isinstance(flines.open_field_lines.solar_feet,
                      coord.SkyCoord)
    assert isinstance(flines.open_field_lines.source_surface_feet,
                      coord.SkyCoord)
    assert isinstance(flines.polarities,
                      np.ndarray)


def test_fline_in_bounds(flines):
    assert np.all(flines[0].coords.radius >= const.R_sun)
    assert np.all(flines[0].coords.radius <= 2.5 * const.R_sun)


def test_rot_warning(dipole_result):
    tracer = tracing.FortranTracer(max_steps=2)
    input, out = dipole_result
    out_frame = out.coordinate_frame
    seed = coord.SkyCoord(0*u.deg, -45*u.deg, 1.01*const.R_sun,
                          frame=out_frame)

    with pytest.warns(UserWarning, match='ran out of steps'):
        tracer.trace(seed, out)

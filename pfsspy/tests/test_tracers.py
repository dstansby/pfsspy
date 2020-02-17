import astropy.coordinates as coord
import astropy.constants as const
import astropy.units as u
import numpy as np
import pytest

import pfsspy
from pfsspy import tracing

from .example_maps import dipole_map, dipole_result


@pytest.mark.parametrize('tracer', [tracing.PythonTracer(),
                                    tracing.FortranTracer()],
                         ids=('python', 'fortran'))
def test_field_lines(dipole_result, tracer):
    input, out = dipole_result
    out_frame = out.coordinate_frame

    seed = coord.SkyCoord(0*u.deg, -45*u.deg, 1.01*const.R_sun, frame=out_frame)
    field_lines = tracer.trace(seed, out)
    assert isinstance(field_lines[0],
                      pfsspy.fieldline.FieldLine)
    assert isinstance(field_lines.open_field_lines.solar_feet,
                      coord.SkyCoord)
    assert isinstance(field_lines.open_field_lines.source_surface_feet,
                      coord.SkyCoord)
    assert isinstance(field_lines.polarities,
                      np.ndarray)


def test_rot_warning(dipole_result):
    tracer = tracing.FortranTracer(max_steps=2)
    input, out = dipole_result
    out_frame = out.coordinate_frame
    seed = coord.SkyCoord(0*u.deg, -45*u.deg, 1.01*const.R_sun,
                          frame=out_frame)

    with pytest.warns(UserWarning, match='ran out of steps'):
        field_lines = tracer.trace(seed, out)

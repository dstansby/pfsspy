import astropy.coordinates as coord
import numpy as np
import pytest

import pfsspy
from pfsspy import tracing

from .example_maps import dipole_map


@pytest.mark.parametrize('tracer', [tracing.PythonTracer(),
                                    tracing.FortranTracer()],
                         ids=('python', 'fortran'))
@pytest.mark.parametrize('seeds', [np.array([1.01, 0, 1.01]),
                                   np.array([[1.01, 0, 1.01],
                                             [1.01, 0, 1.01]])],
                         ids=('1 seed', 'multiple seeds'))
def test_field_lines(dipole_map, seeds, tracer):
    input, out = dipole_map

    field_lines = tracer.trace(seeds, out)
    assert isinstance(field_lines[0],
                      pfsspy.fieldline.FieldLine)
    assert isinstance(field_lines.open_field_lines.solar_feet,
                      coord.SkyCoord)
    assert isinstance(field_lines.open_field_lines.source_surface_feet,
                      coord.SkyCoord)
    assert isinstance(field_lines.polarities,
                      np.ndarray)


def test_rot_warning(dipole_map):
    tracer = tracing.FortranTracer(max_steps=2)
    input, out = dipole_map

    with pytest.warns(UserWarning, match='ran out of steps'):
        field_lines = tracer.trace(np.array([1.01, 0, 1.01]), out)

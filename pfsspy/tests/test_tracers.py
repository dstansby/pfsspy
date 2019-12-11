import astropy.coordinates as coord
import numpy as np
import pytest

import pfsspy
import pfsspy.coords
from pfsspy import tracing


@pytest.fixture
def dipole_map():
    # Test a completely zero input
    ntheta = 30
    nphi = 20
    nr = 10
    rss = 2.5

    phi = np.linspace(0, 2 * np.pi, nphi)
    theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
    theta, phi = np.meshgrid(theta, phi)

    def dipole_Br(r, theta):
        return 2 * np.sin(theta) / r**3

    br = dipole_Br(1, theta).T
    input = pfsspy.Input(br, nr, rss)
    output = pfsspy.pfss(input)
    return input, output


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

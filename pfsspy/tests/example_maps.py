import pytest
import numpy as np
from sunpy.map import Map
from astropy.time import Time
import pfsspy


@pytest.fixture
def zero_map():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))
    header = pfsspy.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    input_map = Map((br, header))

    input = pfsspy.Input(input_map, nr, rss)
    output = pfsspy.pfss(input)
    return input, output


@pytest.fixture
def dipole_map():
    ntheta = 30
    nphi = 20

    phi = np.linspace(0, 2 * np.pi, nphi)
    theta = np.linspace(-np.pi / 2, np.pi / 2, ntheta)
    theta, phi = np.meshgrid(theta, phi)

    def dipole_Br(r, theta):
        return 2 * np.sin(theta) / r**3

    br = dipole_Br(1, theta).T
    header = pfsspy.carr_cea_wcs_header(Time('1992-12-21'), br.shape)
    return Map((br, header))


@pytest.fixture
def dipole_result(dipole_map):
    nr = 10
    rss = 2.5

    input = pfsspy.Input(dipole_map, nr, rss)
    output = pfsspy.pfss(input)
    return input, output

import pytest
import numpy as np
import pfsspy


@pytest.fixture
def zero_map():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))

    input = pfsspy.Input(br, nr, rss)
    output = pfsspy.pfss(input)
    return input, output


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
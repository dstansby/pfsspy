import astropy.units as u
import numpy as np
import pytest

from pfsspy import analytic


@pytest.mark.parametrize('f', [analytic.Bphi, analytic.Btheta])
def test_bangular_rss(f):
    # Check that angular components of the field are zero on source surface
    zss = 2
    f = f(1, 1, zss)
    npoints = 200
    theta = np.linspace(10, 170, npoints) * u.deg
    phi = np.linspace(0, 360, npoints) * u.deg
    assert u.allclose(f(zss, theta, phi), 0)


def test_br_rss():
    zss = 2
    l = 1
    m = 0
    c = zss**(-l-2) * ((2*l + 1) / (l + 1 + l * zss**(-2*l - 1)))
    f = analytic.Br(l, m, zss)
    phi = 0 * u.deg
    theta = 0 * u.deg
    assert f(zss, theta, phi) == 0.5 * np.sqrt(3 / np.pi) * c

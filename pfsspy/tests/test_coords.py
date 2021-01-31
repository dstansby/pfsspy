import astropy.units as u
import numpy as np

import pfsspy.coords


def test_forward_backwards():
    xyz = [1, 0, 0] * u.m
    rho, s, phi = pfsspy.coords.cart2strum(*xyz)
    x1, y1, z1 = pfsspy.coords.strum2cart(rho, s, phi)

    assert u.allclose(xyz[0], x1)
    assert u.allclose(xyz[1], y1)
    assert u.allclose(xyz[2], z1)

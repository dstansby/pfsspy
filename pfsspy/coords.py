r"""
Helper functions for coordinate transformations used in the PFSS domain.
"""

import numpy as np


def sph2cart(r, theta, phi):
    """
    Convert spherical coordinates to cartesian coordinates.
    """
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return x, y, z


def strum2cart(rho, s, phi):
    """
    Convert strumfric coordinates to cartesian coordinates.
    """
    r = np.exp(rho)
    theta = np.arccos(s)
    return sph2cart(r, theta, phi)


def cart2sph(x, y, z):
    """
    Convert cartesian coordinates to spherical coordinates.

    Returns
    -------
    r
    theta
    phi
    """
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    return r, theta, phi


def cart2strum(x, y, z):
    """
    Convert cartesian coordinates to strumfric coordinates.

    Returns
    -------
    rho
    s
    phi
    """
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)
    r = np.sqrt(x**2 + y**2 + z**2)
    s = z / r  # = cos(theta)
    rho = np.log(r)
    return rho, s, phi

r"""
Helper functions for coordinate transformations used in the PFSS domain.

The PFSS solution is calculated on a "strumfric" grid defined by

- :math:`\rho = \log (r)`
- :math:`s = \cos (\theta )`
- :math:`\phi`

where :math:`r, \theta, \phi` are spherical cooridnates that have ranges

- :math:`1 < r < r_{ss}`
- :math:`0 < \theta < \pi`
- :math:`0 < \phi < 2\pi`

The transformation between cartesian coordinates used by the tracer and the
above coordinates is given by

- :math:`x = r\sin (\theta) \cos (\phi)`
- :math:`y = r\sin (\theta) \sin (\phi)`
- :math:`z = r \cos (\theta)`
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


def cart2strum(x, y, z):
    """
    Convert cartesian coordinates to strumfric coordinates.
    """
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)
    r = np.sqrt(x**2 + y**2 + z**2)
    s = z / r  # = cos(theta)
    rho = np.log(r)
    return rho, s, phi

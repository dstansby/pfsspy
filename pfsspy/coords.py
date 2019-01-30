r"""
Helper functions for coordinate transformations used in the PFSS domain.

The PFSS solution is calculated on a "strumfric" grid defined by

- :math:`\rho = \log (r)`
- :math:`s = \cos (\theta )`
- :math:`\phi`

where :math:`r, \theta, \phi` are spherical cooridnates that have ranges

- :math:`1 < r < r_{ss}`
- :math:`-\pi / 2 < \theta < \pi / 2`
- :math:`0 < \phi < 2\pi`

The transformation between cartesian coordinates used by the tracer and the
above coordinates is given by

- :math:`x = r\cos (\theta) \cos (\phi)`
- :math:`y = r\cos (\theta) \sin (\phi)`
- :math:`z = r \sin (\theta)`
"""

import numpy as np


def strum2cart(rho, s, phi):
    """
    Convert strumfric coordinates to cartesian coordinates.
    """
    r = np.exp(rho)
    theta = np.arccos(s)

    x = r * np.cos(theta) * np.cos(phi)
    y = r * np.cos(theta) * np.sin(phi)
    z = r * np.sin(theta)
    return x, y, z


def cart2strum(x, y, z):
    """
    Convert cartesian coordinates to strumfric coordinates.
    """
    phi = (np.arctan2(y, x) + 2 * np.pi) % (2 * np.pi)
    r = np.sqrt(x**2 + y**2 + z**2)
    s = z / r  # = cos(theta)
    rho = np.log(r)
    return rho, s, phi

"""
Helper functions for coordinate transformations used in the PFSS domain.
"""
import astropy.constants as const
from astropy.coordinates import spherical_to_cartesian, cartesian_to_spherical
import astropy.units as u
import numpy as np


__all__ = ['strum2cart', 'cart2strum']


@u.quantity_input
def strum2cart(rho: u.dimensionless_unscaled,
               s: u.dimensionless_unscaled,
               phi: u.rad):
    """
    Convert strumfric coordinates to cartesian coordinates.
    """
    r = np.exp(rho) * const.R_sun
    theta = np.arccos(s)
    return spherical_to_cartesian(r, theta, phi)


@u.quantity_input
def cart2strum(x: u.m, y: u.m, z: u.m):
    """
    Convert cartesian coordinates to strumfric coordinates.

    Returns
    -------
    rho
    s
    phi
    """
    r, theta, phi = cartesian_to_spherical(x, y, z)
    s = np.cos(theta)
    rho = np.log(r / const.R_sun)
    return rho, s, phi

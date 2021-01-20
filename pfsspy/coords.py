"""
Helper functions for coordinate transformations used in the PFSS domain.
"""
import astropy.constants as const
from astropy.coordinates import spherical_to_cartesian, cartesian_to_spherical
import astropy.units as u
import numpy as np


__all__ = ['strum2cart', 'cart2strum']


def strum2cart(rho, s, phi):
    """
    Convert strumfric coordinates to cartesian coordinates.
    """
    r = np.exp(rho) * const.R_sun
    theta = np.arccos(s) * u.rad
    return spherical_to_cartesian(r, theta, phi)


def cart2strum(x, y, z):
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

r"""
Analytic inputs and solutions to the PFSS equations.

This sub-module contains functions to generate solutions to the PFSS equations
in the case where the input field is a single spherical harmonic, specified
with the spherical harmonic numbers ``l, m``.

All angular quantities must be passed as astropy quantities. All radial
quantities are passed normalised to the source surface radius, and therefore
can be passed as normal scalar values.

Angular definitions
-------------------
- ``theta`` is the polar angle, in the range :math:`0, \pi`
  (ie. the co-latitude).
- ``phi`` is the azimuthal angle, in the range :math:`0, 2\pi`.

Using this module requires `sympy` to be installed.
"""

import astropy.units as u
import numpy as np
import scipy.special
import sympy


@u.quantity_input
def _normalise_angles(theta: u.deg, phi: u.deg):
    """
    - Strip units and return radians
    """
    theta = theta.to_value(u.rad)
    phi = phi.to_value(u.rad)
    return theta, phi


def _Ynm(l, m, theta, phi):
    """
    Return values of spherical harmonic with numbers l, m at coordiantes
    theta, phi.
    """
    # Note swapped arguments phi, theta, as scipy has a different
    # definition of these
    return -scipy.special.sph_harm(m, l, phi, theta)


def _cot(theta):
    return 1 / np.tan(theta)


_extras = {'Ynm': _Ynm, 'cot': _cot, 'exp': np.exp}


def _spherical_harmonic_sympy(l, m):
    """
    Return a complex spherical harmonic with numbers ``l, m``.

    Parameters
    ----------
    l, m: int
        Spherical harmonic numbers.

    Returns
    -------
    harm :
    theta, phi : sympy.core.symbol.Symbol

    See also
    --------
    sympy.functions.special.spherical_harmonics.Ynm
    """
    L, M = sympy.symbols('l, m')
    theta, phi = sympy.symbols('theta, phi')
    harm = sympy.Ynm(L, M, theta, phi)
    if m < 0:
        # Phase shift to align definition of Ymn with defnition in paper.
        harm *= -1j
    harm = harm.subs([(L, l), (M, m)])
    return harm, theta, phi


def _c(l, zss):
    """
    """
    def cl(z):
        return (z**(-l - 2) *
                (l + 1 + l * (z / zss)**(2 * l + 1)) /
                (l + 1 + l * zss**(-2 * l - 1)))

    return cl


def _d(l, zss):
    """
    """
    def dl(z):
        return (z**(-l - 2) *
                (1 - (z / zss)**(2 * l + 1)) /
                (l + 1 + l * zss**(-2 * l - 1)))

    return dl


def Br(l, m, zss):
    """
    Analytic radial component of magnetic field on the source surface.

    Parameters
    ----------
    l, m: int
        Spherical harmonic numbers.
    zss: float
        Source surface radius (as a fraction of the solar radius).

    Returns
    -------
    function :
        Has the signature ``Br(z, theta, phi)``.
    """
    sph, t, p = _spherical_harmonic_sympy(l, m)
    sph = sympy.lambdify((t, p), sph, _extras)

    @u.quantity_input
    def f(z, theta: u.deg, phi: u.deg):
        theta, phi = _normalise_angles(theta, phi)
        return _c(l, zss)(z) * np.real(sph(theta, phi))

    return f


def Btheta(l, m, zss):
    """
    Analytic theta component of magnetic field on the source surface.

    Parameters
    ----------
    l, m: int
        Spherical harmonic numbers.
    zss: float
        Source surface radius (as a fraction of the solar radius).

    Returns
    -------
    function :
        Has the signature ``Btheta(z, theta, phi)``.
    """
    sph, t, p = _spherical_harmonic_sympy(l, m)
    sph = sympy.diff(sph, t)
    sph = sympy.lambdify((t, p), sph, [_extras, 'numpy'])

    @u.quantity_input
    def f(z, theta: u.deg, phi: u.deg):
        theta, phi = _normalise_angles(theta, phi)
        return _d(l, zss)(z) * np.real(sph(theta, phi))

    return f


def Bphi(l, m, zss):
    """
    Analytic phi component of magnetic field on the source surface.

    Parameters
    ----------
    l, m: int
        Spherical harmonic numbers.
    zss: float
        Source surface radius (as a fraction of the solar radius).

    Returns
    -------
    function :
        Has the signature ``Bphi(z, theta, phi)``.
    """
    sph, t, p = _spherical_harmonic_sympy(l, m)
    sph = sympy.diff(sph, p)
    sph = sympy.lambdify((t, p), sph, [_extras, 'numpy'])

    @u.quantity_input
    def f(z, theta: u.deg, phi: u.deg):
        theta, phi = _normalise_angles(theta, phi)
        return _d(l, zss)(z) * np.real(sph(theta, phi)) / np.sin(theta)

    return f

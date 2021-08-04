import astropy.units as u
import numpy as np
import sunpy.map
from sympy import acos, asin, cos, lambdify, sin
from sympy.abc import x

import pfsspy.utils
from pfsspy import analytic


def theta_phi(nphi, ns):
    # Return a theta, phi grid with a given numer of points
    phi = np.linspace(0, 2 * np.pi, nphi)
    s = np.linspace(-1, 1, ns)
    s, phi = np.meshgrid(s, phi)
    theta = np.arccos(s)
    return theta * u.rad, phi * u.rad


def pffspy_output(nphi, ns, nrho, rss, l, m):
    # Return the pfsspy solution for given input parameters
    theta, phi = theta_phi(nphi, ns)

    br_in = analytic.Br(l, m, rss)(1, theta, phi)

    header = pfsspy.utils.carr_cea_wcs_header('2020-1-1', br_in.shape)
    input_map = sunpy.map.Map((br_in.T, header))

    pfss_input = pfsspy.Input(input_map, nrho, rss)
    return pfsspy.pfss(pfss_input)


def brss_pfsspy(nphi, ns, nrho, rss, l, m):
    pfsspy_out = pffspy_output(nphi, ns, nrho, rss, l, m)
    return pfsspy_out.bc[0][:, :, -1].T


def brss_analytic(nphi, ns, rss, l, m):
    # Return the analytic solution for given input parameters
    theta, phi = theta_phi(nphi, ns)
    return analytic.Br(l, m, rss)(rss, theta, phi).T


######################
# Field line helpers #
@u.quantity_input
def fr(r: u.m, rss: u.m, l):
    rho = r / rss
    return ((rho**l * (2*l+1)) /
            ((l * rho**(2*l+1)) + l+1))**(1 / (l+1))


flm_dict = {(1, 0): [cos(x), acos(x)],
            (1, 1): [sin(x), asin(x)],
            (2, 1): [cos(2*x)**(1/2), acos(x**2) / 2],
            (2, 2): [sin(x), asin(x)]}

glm_dict = {(1, 1): sin(x) / cos(x),
            (2, 1): (sin(x)**2 / cos(2*x))**1/2}


@u.quantity_input
def theta_fline_coords(r: u.m, rss: u.m, l, m, theta: u.rad):
    """
    r :
        Radial point.
    rss :
        Source surface radius.
    l, m : int
        Spherical harmonic numbers.
    theta :
        Source surface latitude.
    """
    flm = lambdify(x, flm_dict[(l, m)][0], "numpy")
    flm_inv = lambdify(x, flm_dict[(l, m)][1], "numpy")
    theta_out = flm_inv(flm(theta) * fr(r, rss, l))
    theta_out *= np.sign(theta_out) * np.sign(theta)
    return theta_out


@u.quantity_input
def phi_fline_coords(r: u.m, rss: u.m, l, m, theta: u.rad, phi: u.rad):
    """
    r :
        Radial point.
    rss :
        Source surface radius.
    l, m : int
        Spherical harmonic numbers.
    theta, phi :
        Source surface latitude and longitude.
    """
    theta_fline = theta_fline_coords(r, rss, l, m, theta)
    glm = lambdify(x, glm_dict[(l, m)], "numpy")
    phi_out = np.arcsin(glm(theta_fline) *
                        np.sin(phi) /
                        glm(theta))
    pi12 = np.pi / 2 * u.rad
    pi32 = 3 * np.pi / 2 * u.rad
    to_wrap = (pi12 < phi) & (phi < pi32)
    phi_out[to_wrap] = -phi_out[to_wrap] + np.pi * u.rad
    phi_out[pi32 < phi] += 2 * np.pi * u.rad
    return phi_out

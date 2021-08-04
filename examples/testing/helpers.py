import astropy.units as u
import numpy as np
import sunpy.map

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

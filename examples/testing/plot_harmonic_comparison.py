"""
Spherical harmonic comparisons
==============================
Comparing analytical spherical harmonic solutions to PFSS output.
"""
import astropy.units as u
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import numpy as np
import sunpy.map

import pfsspy
from pfsspy import analytic

from helpers import LMAxes


###############################################################################
# Setup some useful functions for testing
def theta_phi(nphi, ns):
    # Return a theta, phi grid with a given numer of points
    phi = np.linspace(0, 2 * np.pi, nphi)
    s = np.linspace(-1, 1, ns)
    s, phi = np.meshgrid(s, phi)
    theta = np.arccos(s)
    return theta * u.rad, phi * u.rad


def brss_pfsspy(nphi, ns, nrho, rss, l, m):
    # Return the pfsspy solution for given input parameters
    theta, phi = theta_phi(nphi, ns)

    br_in = analytic.Br(l, m, rss)(1, theta, phi)

    header = pfsspy.utils.carr_cea_wcs_header('2020-1-1', br_in.shape)
    input_map = sunpy.map.Map((br_in.T, header))

    pfss_input = pfsspy.Input(input_map, nrho, rss)
    pfss_output = pfsspy.pfss(pfss_input)
    return pfss_output.bc[0][:, :, -1].T


def brss_analytic(nphi, ns, rss, l, m):
    # Return the analytic solution for given input parameters
    theta, phi = theta_phi(nphi, ns)
    return analytic.Br(l, m, rss)(rss, theta, phi).T


###############################################################################
# Compare the the pfsspy solution to the analytic solutions. Cuts are taken
# on the source surface at a constant phi value to do a 1D comparison.
nphi = 360
ns = 180
rss = 2
nrho = 20

nl = 5
axs = LMAxes(nl=nl)

for l in range(1, nl+1):
    for m in range(-l, l+1):
        print(f'l={l}, m={m}')
        ax = axs[l, m]

        br_pfsspy = brss_pfsspy(nphi, ns, nrho, rss, l, m)
        br_actual = brss_analytic(nphi, ns, rss, l, m)

        ax.plot(br_pfsspy[:, 15], label='pfsspy')
        ax.plot(br_actual[:, 15], label='analytic')
        if l == 1 and m == 0:
            ax.xaxis.set_major_formatter(mticker.StrMethodFormatter('{x}°'))
            ax.xaxis.set_ticks([0, 90, 180])
            ax.xaxis.tick_top()
            ax.spines['top'].set_visible(True)
        ax.set_xlim(0, 180)
        ax.axhline(0, linestyle='--', linewidth=0.5, color='black')


plt.show()

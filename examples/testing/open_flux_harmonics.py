"""
Open flux
=========
Comparing total unsigned flux to analytic solutions. This is done on a fixed
number of radial grid points, and plotted as a function of spherical harmonic.
"""

###############################################################################
# First, import required modules
import functools
import json

import astropy.units as u
import matplotlib.colors as mcolor
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.integrate
from helpers import brss_pfsspy
from tqdm import tqdm

from pfsspy import analytic

def open_flux_analytic(l, m, zss):
    Br = analytic.Br(l, m, zss)
    Br = functools.partial(Br, zss)
    absBr = lambda theta, phi: np.abs(Br(theta * u.rad, phi *u.rad)) * np.sin(theta)
    res = scipy.integrate.nquad(absBr, ranges=[[0, np.pi], [0, 2 * np.pi]])
    return res[0]


def open_flux_numeric(l, m, zss, nrho):
    nphi = 360
    ns = 180
    br = brss_pfsspy(nphi, ns, nrho, zss, l, m)
    return np.sum(np.abs(br)) * (4 * np.pi) / nphi / ns


###############################################################################
# Set the source surface height, and the (l, m) values to investigate
zss = 2
nrho = 40

results = {}

for l in range(1, 6):
    results[l] = {}
    for m in range(0, l + 1):
        print(f"l={l}, m={m}")
        flux_analytic = open_flux_analytic(l, m, zss)
        flux_numeric = open_flux_numeric(l, m, zss, nrho)
        results[l][m] = flux_numeric / flux_analytic

# open file for writing, "w"
with open("open_flux_harmonics.json", "w") as f:
    # write json object to file
    f.write(json.dumps(results))


with open("open_flux_harmonics.json", "r") as f:
    results = json.load(f, parse_int=int)
print(results)
###############################################################################
# Plot results
fig, ax = plt.subplots()
norm = mcolor.Normalize(vmin=1, vmax=1.05)
for lstr in results:
    l = int(lstr)
    data = np.atleast_2d(list(results[lstr].values())).T
    im = ax.imshow(data, extent=[l-0.5, l+0.5, -0.5, l+0.5],
                   norm=norm)

fig.colorbar(im, extend='both')
ax.set_xlim(0.5, l+0.5)
ax.set_ylim(-0.5, l+0.5)
ax.xaxis.set_major_locator(mticker.MultipleLocator(1))
ax.yaxis.set_major_locator(mticker.MultipleLocator(1))
def fmt(x, pos):
    return str(int(x))
ax.xaxis.set_major_formatter(fmt)
ax.yaxis.set_major_formatter(fmt)
plt.show()

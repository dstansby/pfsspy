"""
Open flux
=========
Comparing total unsigned flux to analytic solutions.

This script calculates both analytic and numerical values, and saves them
to a .json file. This can be read in by ``plot_open_flux_harmonics.py`` to
visualise the result.
"""

###############################################################################
# First, import required modules
import functools
import json
from collections import defaultdict

import astropy.units as u
import numpy as np
import scipy.integrate
from helpers import brss_pfsspy

from pfsspy import analytic


def open_flux_analytic(l, m, zss):
    Br = analytic.Br(l, m, zss)
    Br = functools.partial(Br, zss)

    def absBr(theta, phi):
        return np.abs(Br(theta * u.rad, phi * u.rad)) * np.sin(theta)
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

results = {'numeric': defaultdict(dict),
           'analytic': defaultdict(dict)}

for l in range(1, 6):
    for m in range(-l, l + 1):
        print(f"l={l}, m={m}")
        if -m in results['analytic'][l]:
            # Analytic flux for m = -m is the same
            flux_analytic = results['analytic'][l][-m]
        else:
            flux_analytic = open_flux_analytic(l, m, zss)

        results['analytic'][l][m] = float(flux_analytic)
        flux_numeric = open_flux_numeric(l, m, zss, nrho)
        results['numeric'][l][m] = float(flux_numeric)

# open file for writing, "w"
with open("results/open_flux_harmonics.json", "w") as f:
    # write json object to file
    f.write(json.dumps(results))

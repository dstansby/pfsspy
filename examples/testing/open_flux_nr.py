"""
Open flux and radial grid points (calcuations)
==============================================
Comparing total unsigned flux to analytic solutions.

This script caclulates the ratio of numeric to analytic total unsigned open
fluxes in PFSS solutions of spherical harmonics, as a function of the number of
radial grid cells in the pfsspy grid.
"""
import functools

import astropy.units as u
import numpy as np
import pandas as pd
import scipy.integrate
from helpers import brss_pfsspy
from tqdm import tqdm

from pfsspy import analytic


###############################################################################
# Define a function to calculate the analytic total open flux
def open_flux_analytic(l, m, zss):
    """
    Calculate analytic unsigned open flux for spherical harmonic (*l*, *m*) and
    a source surface radius of *zss*.
    """
    Br = analytic.Br(l, m, zss)
    Br = functools.partial(Br, zss)

    def absBr(theta, phi):
        return np.abs(Br(theta * u.rad, phi * u.rad)) * np.sin(theta)

    res = scipy.integrate.nquad(absBr, ranges=[[0, np.pi], [0, 2 * np.pi]])
    return res[0]


###############################################################################
# Define a function to calculate the numeric total open flux (ie. calculated
# by pfsspy)
def open_flux_numeric(l, m, zss, nrho):
    """
    Calculate numerical unsigned open flux for spherical harmonic (*l*, *m*)
    a source surface radius of *zss* and *nrho* grid points in the radial
    direction.
    """
    nphi = 360
    ns = 180
    br = brss_pfsspy(nphi, ns, nrho, zss, l, m)
    return np.sum(np.abs(br)) * (4 * np.pi) / nphi / ns


###############################################################################
# Set the source surface height and range of radial grid points
zss = 2
nrhos = np.arange(10, 51, 2)
print(f'nrhos={nrhos}')


###############################################################################
# Loop through spherical harmonics and do the calculations. Only the ratio
# of fluxes between the analytic and numeric solutions is saved.
df = pd.DataFrame(index=nrhos, columns=[])

for l in range(1, 6):
    for m in range(-l, l+1):
        lm = str(l) + str(m)
        print(f"l={l}, m={m}")
        flux_analytic = open_flux_analytic(l, m, zss)
        flux_numeric = []
        for nrho in tqdm(nrhos):
            flux_numeric.append(open_flux_numeric(l, m, zss, nrho))
        flux_numeric = np.array(flux_numeric)
        flux_ratio = flux_numeric / flux_analytic
        df[lm] = flux_ratio


###############################################################################
# Save a copy of the data
df.to_csv('results/open_flux_results.csv')

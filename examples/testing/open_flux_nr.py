"""
Open flux
=========
Comparing total unsigned flux to analytic solutions. This is done as a function
of the number of radial grid cells used in the PFSS calculation.
"""

###############################################################################
# First, import required modules
import functools

import astropy.units as u
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
lms = [10]

fig, ax = plt.subplots()
# nrhos = np.arange(10, 51, 2)
nrhos = [10, 50]

df = pd.DataFrame(index=nrhos, columns=lms)

for lm in lms:
    l = lm // 10
    m = lm % 10
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
df.to_csv('open_flux_results.csv')
df = pd.read_csv('open_flux_results.csv', index_col=0)

for lm in lms:
    l = lm // 10
    m = lm % 10
    color = {1: 'tab:blue', 2: 'tab:orange', 3: 'tab:green'}[l]
    marker = {0: 'o', 1: 10, 2: 11}[m]
    ax.plot(nrhos, df[str(lm)],
            marker=marker, color=color,
            label=f'l={l}, m={m}')

ax.legend()
ax.set_ylim(1)
ax.set_xlim(8)
ax.yaxis.grid(linestyle='--')
ax.set_xlabel('$n_{r}$')
ax.set_ylabel('$\Phi_{pfsspy} / \Phi_{analytic}$')
plt.show()

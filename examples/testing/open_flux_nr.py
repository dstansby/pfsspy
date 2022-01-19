"""
Open flux and radial grid points (calcuations)
==============================================
Comparing total unsigned flux to analytic solutions.

This script caclulates the ratio of numeric to analytic total unsigned open
fluxes in PFSS solutions of spherical harmonics, as a function of the number of
radial grid cells in the pfsspy grid.
"""

import numpy as np
import pandas as pd
from tqdm import tqdm

from helpers import open_flux_analytic, open_flux_numeric, result_dir

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
df.to_csv(result_dir / 'open_flux_results.csv')

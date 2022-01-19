"""
Open flux comparison (calculations)
===================================
Comparing total unsigned flux to analytic solutions.

This script calculates both analytic and numerical values of the total unsigned
open flux within PFSS solutions of single spherical harmonics, and saves them
to a .json file. This can be read in by ``plot_open_flux_harmonics.py`` to
visualise the result.
"""
import json
from collections import defaultdict

from helpers import open_flux_analytic, open_flux_numeric

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

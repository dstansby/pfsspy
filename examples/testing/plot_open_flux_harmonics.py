"""
Open flux comparison
====================
Comparing total unsigned flux to analytic solutions.

This script plots results generated by ``open_flux_harmonics.py``, comparing
analytic and numerical values of the total unsigned open flux within PFSS
solutions of single spherical harmonics.
"""
import json

import matplotlib.cm as cm
import matplotlib.colors as mcolor
import matplotlib.pyplot as plt

from helpers import LMAxes

with open("results/open_flux_harmonics.json", "r") as f:
    results = json.load(f, parse_int=int)

axs = LMAxes(nl=5)
norm = mcolor.Normalize(vmin=1, vmax=1.06)
cmap = plt.get_cmap('plasma')

for lstr in results['analytic']:
    for mstr in results['analytic'][lstr]:
        l, m = int(lstr), int(mstr)

        ax = axs[l, m]
        ax.set_facecolor(cmap(norm(results['numeric'][lstr][mstr] /
                                   results['analytic'][lstr][mstr])))
        ax.set_aspect('equal')

cbar = axs.fig.colorbar(cm.ScalarMappable(norm=norm, cmap=cmap),
                        ax=axs.all_axs)
cbar.ax.set_ylabel(r'$\frac{\Phi_{pfsspy}}{\Phi_{analytic}}$', rotation=0,
                   size=18, labelpad=27, va='center')
plt.show()

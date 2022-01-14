"""
Visualising tracer step size
============================
"""

###############################################################################
# First, import required modules
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

import pandas as pd
import numpy as np

from helpers import LMAxes, figdir


nl = 3
axs = LMAxes(nl=nl)

for l in range(1, nl+1):
    for m in range(-l, l+1):
        print(l, m)
        ax = axs[l, m]
        try:
            dphis = pd.read_hdf(f'results/dphis_{l}{m}.hdf', 'table')
            dthetas = pd.read_hdf(f'results/dthetas_{l}{m}.hdf', 'table')
        except FileNotFoundError:
            print(f'❌ Files not found for l={l}, m={m}')
            continue

        for data in [dphis, dthetas]:
            vals = np.ma.masked_array(data.values, np.abs(data.values) > 30)
            ax.plot(data.index, np.nanmax(vals, axis=1), marker='.')

        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.set_ylim(0.5e-1, 2e1)
        for x in [1, 4, 16]:
            ax.axvline(x, color='k', linewidth=1, linestyle='--', alpha=0.2)
        for y in [1e-1, 1, 1e1]:
            ax.axhline(y, color='k', linewidth=1, linestyle='--', alpha=0.2)

        if l == 1 and m == 1:
            ax.xaxis.tick_top()
            ax.yaxis.tick_right()
            ax.xaxis.set_ticks([1, 4, 16])
            ax.xaxis.set_major_formatter(mticker.ScalarFormatter())
            ax.yaxis.set_major_formatter(mticker.StrMethodFormatter('{x}°'))
            ax.xaxis.set_ticks([], minor=True)
            ax.yaxis.set_ticks([], minor=True)
            ax.set_xlabel('Step size')
            ax.set_ylabel('Max\nerror', rotation=0, labelpad=15, va='center')
            ax.xaxis.set_label_position('top')
            ax.yaxis.set_label_position('right')
        else:
            ax.yaxis.set_major_formatter(mticker.NullFormatter())
            for minor in [True, False]:
                ax.xaxis.set_ticks([], minor=minor)
                ax.yaxis.set_ticks([], minor=minor)


axs.fig.savefig(figdir / 'tracer_step_size.pdf', bbox_inches='tight')
plt.show()

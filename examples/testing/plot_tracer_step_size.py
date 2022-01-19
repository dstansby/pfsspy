"""
Tracer step size
================
"""
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import pandas as pd
from helpers import LMAxes

nl = 3

axs = LMAxes(nl=nl)

for l in range(1, nl+1):
    for m in range(-l, l+1):
        ax = axs[l, m]
        try:
            dphis = pd.read_csv(f'results/flines/dphis_{l}{m}.csv',
                                header=None, index_col=0)
            dthetas = pd.read_csv(f'results/flines/dthetas_{l}{m}.csv',
                                  header=None, index_col=0)
            print(l, m)
        except FileNotFoundError:
            print(f'❌ Files not found for l={l}, m={m}')
            continue

        for data in [dphis, dthetas]:
            ax.plot(data.index, data.values, marker='.')

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


plt.show()

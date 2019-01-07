import numpy as np
import matplotlib.pyplot as plt


def radial_cut(phi, costheta, field, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    mesh = ax.pcolormesh(np.rad2deg(phi), costheta, field, cmap='RdBu')
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\cos (\theta)$')
    ax.set_xlim(0, 360)
    ax.set_ylim(-1, 1)

    return mesh

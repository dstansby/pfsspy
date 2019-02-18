import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker


def radial_cut(phi, costheta, field, ax=None):
    if ax is None:
        fig, ax = plt.subplots()

    vlim = np.max(np.abs(field))
    mesh = ax.pcolormesh(np.rad2deg(phi), costheta, field,
                         cmap='RdBu', vmin=-vlim, vmax=vlim)
    ax.set_xlabel(r'$\phi$')
    ax.set_ylabel(r'$\cos (\theta)$')
    ax.set_xlim(0, 360)
    ax.set_ylim(-1, 1)
    ax.set_aspect(0.5 * 360 / 2)

    ax.xaxis.set_major_locator(mticker.MultipleLocator(base=60))
    ax.xaxis.set_minor_locator(mticker.MultipleLocator(base=30))

    return mesh


def contour(phi, costheta, field, levels, ax=None, **kwargs):
    """
    Parameters
    ----------
    phi :
    costheta :
    field :
    levels :
    ax : Axes, optional
        Axes to plot to. If ``None`` a new figure is created.
    **kwargs :
        Keyword arguments are handed to `ax.contour`.
    """
    if ax is None:
        fig, ax = plt.subplots()
    phi, theta = np.meshgrid(phi, costheta)
    ax.contour(np.rad2deg(phi), theta, field, levels=levels, **kwargs)

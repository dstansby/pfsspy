import functools
import numpy as np


class Grid:
    r"""
    Grid on which the pfsspy solution is calculated.

    Notes
    -----
    The PFSS solution is calculated on a "strumfric" grid defined by

    - :math:`\rho = \log (r)`
    - :math:`s = \cos (\theta )`
    - :math:`\phi`

    where :math:`r, \theta, \phi` are spherical cooridnates that have ranges

    - :math:`1 < r < r_{ss}`
    - :math:`0 < \theta < \pi`
    - :math:`0 < \phi < 2\pi`
    """
    def __init__(self, ns, nphi, nr, rss):
        self.ns = ns
        self.nphi = nphi
        self.nr = nr
        self.rss = rss

    @property
    def ds(self):
        """
        Cell size in cos(theta).
        """
        return 2.0 / self.ns

    @property
    def dr(self):
        """
        Cell size in log(r).
        """
        return np.log(self.rss) / self.nr

    @property
    def dp(self):
        """
        Cell size in phi.
        """
        return 2 * np.pi / self.nphi

    @property
    def rc(self):
        """
        Location of the centre of cells in log(r).
        """
        return np.linspace(0.5 * self.dr, np.log(self.rss) - 0.5 * self.dr, self.nr)

    @property
    def sc(self):
        """
        Location of the centre of cells in cos(theta).
        """
        return np.linspace(-1 + 0.5 * self.ds, 1 - 0.5 * self.ds, self.ns)

    @property
    def pc(self):
        """
        Location of the centre of cells in phi.
        """
        return np.linspace(0.5 * self.dp, 2 * np.pi - 0.5 * self.dp, self.nphi)

    @property
    def rg(self):
        """
        Location of the edges of grid cells in log(r).
        """
        return np.linspace(0, np.log(self.rss), self.nr + 1)

    @property
    def sg(self):
        """
        Location of the edges of grid cells in cos(theta).
        """
        return np.linspace(-1, 1, self.ns + 1)

    @property
    def pg(self):
        """
        Location of the edges of grid cells in phi.
        """
        return np.linspace(0, 2 * np.pi, self.nphi + 1)

    @property
    def _grid_spacing(self):
        """
        Return grid spacing as a 3-len list.
        """
        return [self.dp, self.ds, self.dr]

    @property
    @functools.lru_cache()
    def _sqrtsg_correction(self):
        """
        The sqrt(1 - sg**2) correction needed to trace natively. Computed here
        once and cached for performance.
        """
        # Correct s direction for coordinate system distortion
        _, sg, _ = np.meshgrid(self.pg, self.sg, self.rg,
                               indexing='ij')
        return np.sqrt(1 - sg**2)

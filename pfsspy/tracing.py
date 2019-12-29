import abc

import astropy.constants as const
import numpy as np
import sunpy.coordinates.frames as frames

import pfsspy
import pfsspy.fieldline as fieldline


class Tracer(abc.ABC):
    """
    Abstract base class for a streamline tracer.
    """
    @abc.abstractmethod
    def trace(self, seeds, output):
        """
        Parameters
        ----------
        seeds : (n, 3) array
            Coordinaes of the magnetic field seed points, in cartesian
            coordinates.
        output : pfsspy.Output
            pfss output.

        Returns
        -------
        streamlines : FieldLines
            Traced field lines.
        """
        pass

    @staticmethod
    def validate_seeds_shape(seeds):
        """
        Check that *seeds* has the right shape.
        """
        if not seeds.ndim == 2:
            raise ValueError(f'seeds must be a 2D array (got shape {seeds.shape})')
        if not seeds.shape[1] == 3:
            raise ValueError(f'seeds must be a (n, 3) shaped array (got shape {seeds.shape})')

    @staticmethod
    def cartesian_to_coordinate():
        """
        Convert cartesian coordinate outputted by a tracer to a `FieldLine`
        object.
        """


class FortranTracer(Tracer):
    r"""
    Tracer using Fortran code.

    Parameters
    ----------
    max_steps: int, optional
        Maximum number of steps each streamline can take before stopping.
    step_size : float, optional
        Step size as a fraction of cell size at the equator.

    Notes
    -----
    Because the stream tracing is done in spherical coordinates, there is a
    singularity at the poles (ie. :math:`s = \pm 1`), which means seeds placed
    directly on the poles will not go anywhere.
    """
    def __init__(self, max_steps=1000, step_size=0.05):
        from streamtracer import StreamTracer
        self.max_steps = max_steps
        self.step_size = step_size
        self.tracer = StreamTracer(max_steps, step_size)

    def trace(self, seeds, output):
        from streamtracer import VectorGrid
        seeds = np.atleast_2d(seeds)
        self.validate_seeds_shape(seeds)

        # The indexing order on the last index is (phi, s, r)
        vectors = output.bg

        # Correct s direction for coordinate system distortion
        _, sg, _ = np.meshgrid(output.grid.pg, output.grid.sg, output.grid.rg,
                               indexing='ij')
        sqrtsg = np.sqrt(1 - sg**2)
        # phi correction
        with np.errstate(invalid='ignore'):
            vectors[..., 0] /= sqrtsg
        # Technically where s=0 Bphi is now infinite, but because this is
        # singular and Bphi doesn't matter, just set it to a large number
        vectors[sg[..., 0] == 0, 0] = 0.8e+308
        # s correction
        vectors[..., 1] *= -sqrtsg

        grid_spacing = output.grid._grid_spacing
        # Cyclic only in the phi direction
        # (theta direction becomes singular at the poles so it is not cyclic)
        cyclic = [True, False, False]
        origin_coord = [0, -1, 0]
        vector_grid = VectorGrid(vectors, grid_spacing, cyclic=cyclic,
                                 origin_coord=origin_coord)

        # Transform seeds from cartesian to strum
        seeds = pfsspy.coords.cart2strum(seeds[:, 0], seeds[:, 1], seeds[:, 2])
        seeds = np.stack(seeds, axis=-1)
        # Put the seeds in (phi, s, rho) order
        seeds = seeds[:, ::-1]

        # Do the tracing
        self.tracer.trace(seeds, vector_grid)
        xs = self.tracer.xs

        xs = [np.stack(pfsspy.coords.strum2cart(x[:, 2], x[:, 1], x[:, 0]), axis=-1) for x in xs]
        flines = [fieldline.FieldLine(x[:, 0], x[:, 1], x[:, 2], output.dtime, output) for x in xs]
        return fieldline.FieldLines(flines)


class PythonTracer(Tracer):
    """
    Tracer using native python code.

    Uses `scipy.integrate.solve_ivp`, with an LSODA method.
    """
    def __init__(self, atol=1e-4, rtol=1e-4):
        """
        dtf : float, optional
            Absolute tolerance of the tracing.
        rtol : float, optional
            Relative tolerance of the tracing.
        """
        self.atol = atol
        self.rtol = rtol

    def trace(self, seeds, output):
        seeds = np.atleast_2d(seeds)
        self.validate_seeds_shape(seeds)
        flines = []
        for seed in seeds:
            xforw = output._integrate_one_way(1, seed, self.rtol, self.atol)
            xback = output._integrate_one_way(-1, seed, self.rtol, self.atol)
            xback = np.flip(xback, axis=1)
            xout = np.row_stack((xback.T, xforw.T))
            fline = fieldline.FieldLine(xout[:, 0],
                                        xout[:, 1],
                                        xout[:, 2],
                                        output.dtime,
                                        output)
            flines.append(fline)
        return fieldline.FieldLines(flines)

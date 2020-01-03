import abc
import warnings

import astropy.constants as const
import astropy.coordinates as astrocoords
import astropy.units as u
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
    def validate_seeds(seeds):
        """
        Check that *seeds* has the right shape and is the correct type.
        """
        if not isinstance(seeds, astrocoords.SkyCoord):
            raise ValueError('seeds must be SkyCoord object '
                             f'(got {type(seeds)} instead)')

    @staticmethod
    def cartesian_to_coordinate():
        """
        Convert cartesian coordinate outputted by a tracer to a `FieldLine`
        object.
        """

    @staticmethod
    def transform_seeds(seeds, output):
        """
        Transform *seeds* to the coordinate system of *output*.

        Parameters
        ----------
        seeds : astropy.coordinates.SkyCoord
        output : pfsspy.Output
        """
        return seeds.transform_to(output.coordinate_frame)


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
    def __init__(self, max_steps=1000, step_size=0.01):
        try:
            from streamtracer import StreamTracer
        except ModuleNotFoundError as e:
            raise RuntimeError(
                'Using FortranTracer requires the streamtracer module, '
                'but streamtracer could not be loaded '
                f'with the following error:\n{e}')
        self.max_steps = max_steps
        self.step_size = step_size
        self.tracer = StreamTracer(max_steps, step_size)

    @staticmethod
    def vector_grid(output):
        """
        Create a `streamtracer.VectorGrid` object from an `~pfsspy.Ouput`.
        """
        from streamtracer import VectorGrid

        # The indexing order on the last index is (phi, s, r)
        vectors = output.bg.copy()

        # Correct s direction for coordinate system distortion
        sqrtsg = output.grid._sqrtsg_correction
        # phi correction
        with np.errstate(divide='ignore', invalid='ignore'):
            vectors[..., 0] /= sqrtsg
        # Technically where s=0 Bphi is now infinite, but because this is
        # singular and Bphi doesn't matter, just set it to a large number
        vectors[~np.isfinite(vectors[..., 0]), 0] = 0.8e+308
        # s correction
        vectors[..., 1] *= -sqrtsg

        grid_spacing = output.grid._grid_spacing
        # Cyclic only in the phi direction
        # (theta direction becomes singular at the poles so it is not cyclic)
        cyclic = [True, False, False]
        origin_coord = [0, -1, 0]
        vector_grid = VectorGrid(vectors, grid_spacing, cyclic=cyclic,
                                 origin_coord=origin_coord)
        return vector_grid

    def trace(self, seeds, output):
        self.validate_seeds(seeds)
        seeds = self.transform_seeds(seeds, output)

        phi = seeds.lon + 180 * u.deg
        phi.wrap_angle = 360 * u.deg
        phi = phi.to_value(u.rad)
        if not np.all((0 <= phi) & (phi <= 2 * np.pi)):
            raise ValueError('Some phi coords not in range [0, 2pi]')

        s = np.sin(seeds.lat).to_value(u.dimensionless_unscaled)
        if not np.all((-1 <= s) & (s <= 1)):
            raise ValueError('Some s coords not in range [-1, 1]')

        rho = np.log(seeds.radius.to_value(const.R_sun))

        seeds = np.atleast_2d(np.stack((phi, s, rho), axis=-1))

        # Get a grid
        vector_grid = self.vector_grid(output)

        # Do the tracing
        self.tracer.trace(seeds, vector_grid)
        xs = self.tracer.xs
        rots = self.tracer.ROT
        if np.any(rots == 1):
            warnings.warn(
                'At least one field line ran out of steps during tracing.\n'
                'You should probably increase max_steps '
                f'(currently set to {self.max_steps}) and try again.')

        xs = [np.stack(pfsspy.coords.strum2cart(x[:, 2], x[:, 1], x[:, 0]), axis=-1) for x in xs]
        # Hacky way to rotate back by 180deg
        for xout in xs:
            xout[:, 0:2] *= -1
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
        self.validate_seeds(seeds)
        seeds = self.transform_seeds(seeds, output)
        seeds.representation_type = 'cartesian'
        x = -seeds.x.to_value(const.R_sun)
        y = -seeds.y.to_value(const.R_sun)
        z = seeds.z.to_value(const.R_sun)

        seeds = np.atleast_2d(np.stack((x, y, z), axis=-1))

        flines = []
        for seed in seeds:
            xforw = output._integrate_one_way(1, seed, self.rtol, self.atol)
            xback = output._integrate_one_way(-1, seed, self.rtol, self.atol)
            xback = np.flip(xback, axis=1)
            xout = np.row_stack((xback.T, xforw.T))
            # Hacky way to roate back by 180deg
            xout[:, 0:2] *= -1
            fline = fieldline.FieldLine(xout[:, 0],
                                        xout[:, 1],
                                        xout[:, 2],
                                        output.dtime,
                                        output)
            flines.append(fline)
        return fieldline.FieldLines(flines)

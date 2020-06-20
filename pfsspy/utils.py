import sunpy.map
import sunpy.io
import numpy as np
import os
import warnings
import copy
from scipy import interpolate
from astropy import units as u
import astropy.time
import astropy.coordinates as coord


def load_adapt(adapt_path):
    """
    Parse adapt .fts file as a `sunpy.map.MapSequence`

    ADAPT magnetograms contain 12 realizations and their data
    attribute consists of a 3D data cube where each slice is
    the data corresponding to a separate realization of the
    magnetogram. This function loads the raw fits file and
    parses it to a `sunpy.map.MapSequence` object containing
    a `sunpy.map.Map` instance for each realization.

    Parameters
    ----------
    adapt_path : `str`
        Filepath corresponding to an ADAPT .fts file

    Returns
    -------
    adaptMapSequence : `sunpy.map.MapSequence`
    """
    adapt_fits = sunpy.io.fits.read(adapt_path)
    assert adapt_fits[0].header.get('MODEL') == 'ADAPT', \
        f"{os.path.basename(adapt_path)} header['MODEL'] is not 'ADAPT' "
    data_header_pairs = [(map_slice, adapt_fits[0].header)
                         for map_slice in adapt_fits[0].data]
    adaptMapSequence = sunpy.map.Map(data_header_pairs, sequence=True)
    return adaptMapSequence


def carr_cea_wcs_header(dtime, shape):
    """
    Create a Carrington WCS header for a Cylindrical Equal Area (CEA)
    projection. See [1]_ for information on how this is constructed.

    dtime : datetime, None
        Datetime to associate with the map.
    shape : tuple
        Map shape. The first entry should be number of points in longitude, the
        second in latitude.

    References
    ----------
    .. [1] W. T. Thompson, "Coordinate systems for solar image data",
       https://doi.org/10.1051/0004-6361:20054262
    """
    # If datetime is None, put in a dummy value here to make
    # make_fitswcs_header happy, then strip it out at the end
    obstime = dtime or astropy.time.Time('2000-1-1')

    frame_out = coord.SkyCoord(
        0 * u.deg, 0 * u.deg, obstime=obstime,
        frame="heliographic_carrington", observer='sun')
    # Construct header
    header = sunpy.map.make_fitswcs_header(
        shape, frame_out,
        scale=[360 / shape[0],
               180 / shape[1]] * u.deg / u.pix,
        reference_pixel=[(shape[0] / 2) + 0.5, (shape[1] / 2) + 0.5] * u.pix,
        projection_code="CEA")

    # Fill in these missing values
    header['PV1_1'] = 1
    header['PV2_1'] = 1
    # Fix CELT for lat axis
    header['CDELT2'] = (180 / np.pi) * (2 / shape[1])
    # pop out the time if it isn't supplied
    if dtime is None:
        header.pop('date-obs')
    return header


def is_cea_map(m, error=False):
    """
    Returns `True` if *m* is in a cylindrical-equal-area projeciton.

    Parameters
    ----------
    error : bool
        If `True`, raise an error if *m* is not a CEA projection
    """
    for i in ('1', '2'):
        proj = m.meta[f'ctype{i}'][5:8]
        if proj != 'CEA':
            if error:
                raise ValueError(f'Projection type in CTYPE{i} keyword '
                                 f'must be CEA (got "{proj}")')
            return False
    return True


def is_full_sun_synoptic_map(m, error=False):
    """
    Returns `True` if *m* is a synoptic map spanning the solar surface.

    Parameters
    ----------
    error : bool
        If `True`, raise an error if *m* does not span the whole solar surface.
    """
    shape = m.data.shape

    dphi = m.meta['cdelt1']
    phi = shape[1] * dphi
    if not np.allclose(phi, 360, atol=0.1):
        if error:
            raise ValueError('Number of points in phi direction times '
                             'CDELT1 must be close to 360 degrees. '
                             f'Instead got {dphi} x {shape[0]} = {phi}')
        return False

    dtheta = m.meta['cdelt2']
    theta = shape[0] * dtheta * np.pi / 2
    if not np.allclose(theta, 180, atol=0.1):
        if error:
            raise ValueError('Number of points in theta direction times '
                             'CDELT2 times pi/2 must be close to '
                             '180 degrees. '
                             f'Instead got {dtheta} x {shape[0]} * pi/2 = {theta}')
        return False

    return True

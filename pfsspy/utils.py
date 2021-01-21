import os

import numpy as np

import astropy.coordinates as coord
import astropy.time
from astropy import units as u

import sunpy.io
import sunpy.map
import sunpy.time

import pfsspy.map


def fix_hmi_meta(hmi_map):
    """
    Fix non-compliant FITS metadata in HMI maps.

    This function:
        - Corrects CUNIT2 from 'Sine Latitude' to 'deg'
        - Corrects CDELT1 and CDELT2 for a CEA projection
        - Populates the DATE-OBS keyword from the T_OBS keyword
        - Sets the observer coordinate to the Earth

    Notes
    -----
    If you have sunpy > 2.1 installed, this function is not needed as sunpy
    will automatically make these fixes.
    """
    if not isinstance(hmi_map, sunpy.map.sources.HMIMap):
        raise ValueError('Input must be of type HMIMap. '
                         'n.b. if you have sunpy 2.1 installed, '
                         'this function is redundant '
                         'as sunpy 2.1 automatically fixes HMI metadata.')

    if hmi_map.meta['cunit1'] == 'Degree':
        hmi_map.meta['cunit1'] = 'deg'

    if hmi_map.meta['cunit2'] == 'Sine Latitude':
        hmi_map.meta['cunit2'] = 'deg'

        # Since, this map uses the cylindrical equal-area (CEA) projection,
        # the spacing should be modified to 180/pi times the original value
        # Reference: Section 5.5, Thompson 2006
        hmi_map.meta['cdelt2'] = 180 / np.pi * hmi_map.meta['cdelt2']
        hmi_map.meta['cdelt1'] = np.abs(hmi_map.meta['cdelt1'])

    if 'date-obs' not in hmi_map.meta and 't_obs' in hmi_map.meta:
        hmi_map.meta['date-obs'] = sunpy.time.parse_time(hmi_map.meta['t_obs']).isot

    # Fix observer coordinate
    if 'hglt_obs' not in hmi_map.meta:
        hmi_map.meta.update(pfsspy.map._earth_obs_coord_meta(hmi_map.meta['date-obs']))


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
    header = adapt_fits[0].header
    if header['MODEL'] != 'ADAPT':
        raise ValueError(f"{os.path.basename(adapt_path)} header['MODEL'] "
                         "is not 'ADAPT'.")

    data_header_pairs = [(map_slice, header)
                         for map_slice in adapt_fits[0].data]
    adaptMapSequence = sunpy.map.Map(data_header_pairs, sequence=True)
    return adaptMapSequence


def carr_cea_wcs_header(dtime, shape):
    """
    Create a Carrington WCS header for a Cylindrical Equal Area (CEA)
    projection. See [1]_ for information on how this is constructed.

    Parameters
    ----------
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


def _get_projection(m, i):
    return m.meta[f'ctype{i}'][5:8]


def _check_projection(m, proj_code, error=False):
    for i in ('1', '2'):
        proj = _get_projection(m, i)
        if proj != proj_code:
            if error:
                raise ValueError(f'Projection type in CTYPE{i} keyword '
                                 f'must be {proj_code} (got "{proj}")')
            return False
    return True


def is_cea_map(m, error=False):
    """
    Returns `True` if *m* is in a cylindrical equal area projeciton.

    Parameters
    ----------
    m : sunpy.map.GenericMap
    error : bool
        If `True`, raise an error if *m* is not a CEA projection.
    """
    return _check_projection(m, 'CEA', error=error)


def is_car_map(m, error=False):
    """
    Returns `True` if *m* is in a plate carée projeciton.

    Parameters
    ----------
    m : sunpy.map.GenericMap
    error : bool
        If `True`, raise an error if *m* is not a CAR projection.
    """
    return _check_projection(m, 'CAR', error=error)


def is_full_sun_synoptic_map(m, error=False):
    """
    Returns `True` if *m* is a synoptic map spanning the solar surface.

    Parameters
    ----------
    m : sunpy.map.GenericMap
    error : bool
        If `True`, raise an error if *m* does not span the whole solar surface.
    """
    projection = _get_projection(m, 1)
    checks = {'CEA': _is_full_sun_cea,
              'CAR': _is_full_sun_car}
    if projection not in checks.keys():
        raise NotImplementedError('is_full_sun_synoptic_map is only '
                                  'implemented for '
                                  f'{[key for key in checks.keys()]} '
                                  'projections.')
    return checks[projection](m, error)


def _is_full_sun_car(m, error=False):
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
    theta = shape[0] * dtheta
    if not np.allclose(theta, 180, atol=0.1):
        if error:
            raise ValueError('Number of points in theta direction times '
                             'CDELT2 must be close to 180 degrees. '
                             f'Instead got {dtheta} x {shape[0]} = {theta}')
        return False
    return True


def _is_full_sun_cea(m, error=False):
    shape = m.data.shape

    dphi = m.meta['cdelt1']
    phi = shape[1] * dphi
    if not np.allclose(phi, 360, atol=0.1):
        if error:
            raise ValueError('Number of points in phi direction times '
                             'CDELT1 must be close to 360 degrees. '
                             f'Instead got {dphi} x {shape[1]} = {phi}')
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


def car_to_cea(m, method='interp'):
    """
    Reproject a plate-carée map in to a cylindrical-equal-area map.

    The solver used in pfsspy requires a magnetic field map with values
    equally spaced in sin(lat) (ie. a CEA projection), but some maps are
    provided equally spaced in lat (ie. a CAR projection). This function
    reprojects a CAR map into a CEA map so it can be used with pfsspy.

    Parameters
    ----------
    m : sunpy.map.GenericMap
        Input map
    method : str
        Reprojection method to use. Can be ``'interp'`` (default),
        ``'exact'``, or ``'adaptive'``. See :mod:`reproject` for a
        description of the different methods. Note that different methods will
        give different results, and not all will conserve flux.

    Returns
    -------
    output_map : sunpy.map.GenericMap
        Re-projected map. All metadata is preserved, apart from CTYPE{1,2} and
        CDELT2 which are updated to account for the new projection.

    See also
    --------
    :mod:`reproject` for the methods that perform the reprojection.
    """
    from astropy.wcs import WCS
    from reproject import reproject_interp, reproject_exact, reproject_adaptive
    methods = {'adaptive': reproject_adaptive,
               'interp': reproject_interp,
               'exact': reproject_exact}
    if method not in methods:
        raise ValueError(f'method must be one of {methods.keys()} '
                         f'(got {method})')
    reproject = methods[method]
    # Check input map is valid
    is_full_sun_synoptic_map(m, error=True)
    is_car_map(m, error=True)

    header_out = m.wcs.to_header()
    header_out['CTYPE1'] = header_out['CTYPE1'][:5] + 'CEA'
    header_out['CTYPE2'] = header_out['CTYPE2'][:5] + 'CEA'
    header_out['CDELT2'] = 180 / np.pi * 2 / m.data.shape[0]
    wcs_out = WCS(header_out, fix=False)
    # Check if we need to add the heliographic observer (in sunpy <2.1)
    if hasattr(m.wcs, 'heliographic_observer'):
        wcs_out.heliographic_observer = m.observer_coordinate
    data_out = reproject(m, wcs_out, shape_out=m.data.shape,
                         return_footprint=False)

    meta_out = m.meta.copy()
    for key in ['CTYPE1', 'CTYPE2', 'CDELT2']:
        meta_out[key] = header_out[key]
    m_out = sunpy.map.Map(data_out, header_out)
    return m_out

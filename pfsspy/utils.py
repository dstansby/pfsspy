import os

import astropy.constants as const
import astropy.coordinates as coord
import astropy.time
import numpy as np
import sunpy.io
import sunpy.map
import sunpy.time
from astropy import units as u
from astropy.wcs import WCS

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
        hmi_map.meta['date-obs'] = sunpy.time.parse_time(
            hmi_map.meta['t_obs']).isot

    # Fix observer coordinate
    if 'hglt_obs' not in hmi_map.meta:
        hmi_map.meta.update(pfsspy.map._earth_obs_coord_meta(
            hmi_map.meta['date-obs']))


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
        0 * u.deg, 0 * u.deg, const.R_sun, obstime=obstime,
        frame="heliographic_carrington", observer='self')
    # Construct header
    header = sunpy.map.make_fitswcs_header(
        [shape[1], shape[0]], frame_out,
        scale=[360 / shape[0],
               180 / shape[1]] * u.deg / u.pix,
        reference_pixel=[(shape[0] / 2) - 0.5,
                         (shape[1] / 2) - 0.5] * u.pix,
        projection_code="CEA")

    # Fix CDELT for lat axis
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

    dphi = m.scale.axis1
    phi = shape[1] * u.pix * dphi
    if not np.allclose(np.abs(phi), 360 * u.deg, atol=0.1 * u.deg):
        if error:
            raise ValueError('Number of points in phi direction times '
                             'CDELT1 must be close to 360 degrees. '
                             f'Instead got {dphi} x {shape[0]} = {phi}')
        return False

    dtheta = m.scale.axis2
    theta = shape[0] * u.pix * dtheta
    if not np.allclose(theta, 180 * u.deg, atol=0.1 * u.deg):
        if error:
            raise ValueError('Number of points in theta direction times '
                             'CDELT2 must be close to 180 degrees. '
                             f'Instead got {dtheta} x {shape[0]} = {theta}')
        return False
    return True


def _is_full_sun_cea(m, error=False):
    shape = m.data.shape

    dphi = m.scale.axis1
    phi = shape[1] * u.pix * dphi
    if not np.allclose(np.abs(phi), 360 * u.deg, atol=0.1 * u.deg):
        if error:
            raise ValueError('Number of points in phi direction times '
                             'CDELT1 must be close to 360 degrees. '
                             f'Instead got {dphi} x {shape[1]} = {phi}')
        return False

    dtheta = m.scale.axis2
    theta = shape[0] * u.pix * dtheta * np.pi / 2
    if not np.allclose(theta, 180 * u.deg, atol=0.1 * u.deg):
        if error:
            raise ValueError('Number of points in theta direction times '
                             'CDELT2 times pi/2 must be close to '
                             '180 degrees. '
                             f'Instead got {dtheta} x {shape[0]} * pi/2 '
                             f'= {theta}')
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
    # Add reproject import here to avoid import dependency
    from reproject import reproject_adaptive, reproject_exact, reproject_interp
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

    # Create output FITS header
    header_out = m.wcs.to_header()
    header_out['CTYPE1'] = header_out['CTYPE1'][:5] + 'CEA'
    header_out['CTYPE2'] = header_out['CTYPE2'][:5] + 'CEA'
    header_out['CDELT2'] = 180 / np.pi * 2 / m.data.shape[0]
    wcs_out = WCS(header_out, fix=False)

    # Reproject data
    data_out = reproject(m, wcs_out, shape_out=m.data.shape,
                         return_footprint=False)

    return sunpy.map.Map(data_out, header_out)


@u.quantity_input
def roll_map(m, lh_edge_lon: u.deg = 0.0 * u.deg, method='interp'):
    """
    Roll an input synoptic map so that it's left edge corresponds to a specific
    Carrington longitude.

    Roll is performed by changing the FITS header parameter "CRVAL1"
    to the new value of the reference pixel (FITS header parameter CRPIX1)
    corresponding to aligning the left hand edge of the map with
    lh_edge_lon. The altered header is provided to reproject to produce
    the new map.

    Parameters
    ----------
    m : sunpy.map.GenericMap
        Input map
    lh_edge_lon : float
        Desired Carrington longitude (degrees) for left hand edge of map.
        Default is 0.0 which results in a map with the edges at 0/360 degrees
        Carrington  longitude. Input value must be in the range [0,360]
    method : str
        Reprojection method to use. Can be ``'interp'`` (default),
        ``'exact'``, or ``'adaptive'``. See :mod:`reproject` for a
        description of the different methods. Note that different methods will
        give different results, and not all will conserve flux.

    Returns
    -------
    output_map : sunpy.map.GenericMap
        Re-projected map. All metadata is preserved, apart from CRVAL1  which
        encodes the longitude of the reference pixel in the image, and which
        is updated to produce the correct roll.

    See also
    --------
    :mod:`reproject` for the methods that perform the reprojection.
    """
    # Add reproject import here to avoid import dependency
    from reproject import reproject_adaptive, reproject_exact, reproject_interp
    methods = {'adaptive': reproject_adaptive,
               'interp': reproject_interp,
               'exact': reproject_exact}
    if method not in methods:
        raise ValueError(f'method must be one of {methods.keys()} '
                         f'(got {method})')
    if lh_edge_lon > 360.0 * u.deg or lh_edge_lon < 0.0 * u.deg:
        raise ValueError(f"lh_edge_lon must be in the range [0,360])")

    reproject = methods[method]
    # Check input map is valid
    is_full_sun_synoptic_map(m, error=True)

    # Create output FITS header
    header_out = m.wcs.to_header()
    # Note half pixel shift to ensure LH edge leftmost pixel is
    # correctly aligned with map LH edge
    header_out['CRVAL1'] = (lh_edge_lon - 0.5 * header_out['CDELT1'] * u.deg +
                            header_out['CRPIX1'] * header_out['CDELT1'] *
                            u.deg).to_value(u.deg) % 360
    wcs_out = WCS(header_out, fix=False)

    # Reproject data
    data_out = reproject(m, wcs_out, shape_out=m.data.shape,
                         return_footprint=False)

    return sunpy.map.Map(data_out, header_out)

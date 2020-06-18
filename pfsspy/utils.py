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

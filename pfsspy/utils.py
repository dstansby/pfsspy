import sunpy.map
import sunpy.io
import numpy as np
import os
import warnings
import copy
from scipy import interpolate
from astropy import units as u


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

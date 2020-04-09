"""
GONG helper functions
=====================
"""

import os

import numpy as np


def get_gong_map():
    """
    Automatically download and unzip a sample GONG synoptic map.
    """
    if not os.path.exists('190310t0014gong.fits') and not os.path.exists('190310t0014gong.fits.gz'):
        import urllib.request
        urllib.request.urlretrieve(
            'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/mrzqs190310t0014c2215_333.fits.gz',
            '190310t0014gong.fits.gz')

    if not os.path.exists('190310t0014gong.fits'):
        import gzip
        with gzip.open('190310t0014gong.fits.gz', 'rb') as f:
            with open('190310t0014gong.fits', 'wb') as g:
                g.write(f.read())

    return '190310t0014gong.fits'


def fix_gong_header(header):
    """
    Fix various issues with GONG FITS metadata.
    """
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    # Instead of the spacing in sin(lat), this should be 180/pi times that
    # value (see Thompson 2005)
    header['CDELT2'] = 180 / np.pi * header['CDELT2']
    # sunpy complains if the observer is not set; this doesn't significantly
    # affect anything for a heliogrpahic map, so just set the observer at the
    # center of the Sun.
    header['DSUN_OBS'] = 0
    header['HGLN_OBS'] = 0
    header['HGLT_OBS'] = 0
    return header

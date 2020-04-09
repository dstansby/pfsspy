"""
GONG helper functions
=====================
"""

import os

import numpy as np


def get_gong_map():
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
    header['CUNIT1'] = 'deg'
    header['CUNIT2'] = 'deg'
    header['CDELT2'] = 180 / np.pi * header['CDELT2']
    header['DSUN_OBS'] = 0
    header['HGLN_OBS'] = 0
    header['HGLT_OBS'] = 0
    return header

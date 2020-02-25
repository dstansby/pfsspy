import gzip
import os
import urllib.request


def get_gong_map():
    fname = '190310t0014gong.fits'
    if not os.path.exists(fname) and not os.path.exists(f'{fname}.gz'):
        urllib.request.urlretrieve(
            'https://gong2.nso.edu/oQR/zqs/201903/mrzqs190310/mrzqs190310t0014c2215_333.fits.gz',
            f'{fname}.gz')

    if not os.path.exists(fname):
        with gzip.open(f'{fname}.gz', 'rb') as f:
            with open(fname, 'wb') as g:
                g.write(f.read())

    return fname

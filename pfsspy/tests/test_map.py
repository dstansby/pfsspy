import os

import sunpy.map
import pytest

import pfsspy.map


@pytest.fixture
def gong_map():
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


def test_gong_source(gong_map):
    m = sunpy.map.Map(gong_map)
    assert isinstance(m, pfsspy.map.GongSynopticMap)

def test_adapt_map():
    import sunpy.io
    adapt_fits = sunpy.io.fits.read("ftp://gong.nso.edu/adapt/maps/gong/2020/adapt40311_03k012_202001010000_i00005600n1.fts.gz")
    for map_slice in adapt_fits[0].data :
        m = sunpy.map.Map((map_slice,adapt_fits[0].header))
        assert isinstance(m, pfsspy.map.ADAPTMap)

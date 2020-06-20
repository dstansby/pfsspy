import os
import numpy as np

import sunpy.map
import pytest

import pfsspy
from pfsspy import utils
from .example_maps import dipole_map


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


@pytest.fixture
def adapt_map():
    import urllib.request
    urllib.request.urlretrieve(
        "https://gong.nso.edu/adapt/maps/gong/2020/adapt40311_03k012_202001010000_i00005600n1.fts.gz",
        "adapt20200101.fts.gz")

    if not os.path.exists("adapt20200101.fts"):
        import gzip
        with gzip.open('adapt20200101.fts.gz', 'rb') as f:
            with open('adapt20200101.fts', 'wb') as g:
                g.write(f.read())

    return 'adapt20200101.fts'


def test_load_adapt(adapt_map):
    adaptMapSequence = utils.load_adapt(adapt_map)
    assert isinstance(adaptMapSequence, sunpy.map.MapSequence)
    for map_ in adaptMapSequence:
        assert map_.meta['model'] == "ADAPT"


def test_header_generation():
    dtime = None
    ntheta = 180
    nphi = 360
    shape = [nphi, ntheta]
    header = pfsspy.utils.carr_cea_wcs_header(dtime, shape)
    assert header['LONPOLE'] == 0
    assert header['CTYPE1'] == 'CRLN-CEA'
    assert header['CTYPE2'] == 'CRLT-CEA'
    assert header['PV1_1'] == 1
    assert header['PV2_1'] == 1
    assert header['CDELT1'] == 360 / nphi
    np.testing.assert_almost_equal(
        header['CDELT2'], (180 / np.pi) * (2 / ntheta))

    # Extra + 1s are for FITS counting from 1 indexing
    assert header['CRPIX1'] == (nphi / 2) + 0.5 + 1
    assert header['CRPIX2'] == (ntheta / 2) + 0.5 + 1
    assert header['CRVAL1'] == 0
    assert header['CRVAL2'] == 0
    assert header['CUNIT1'] == 'deg'
    assert header['CUNIT2'] == 'deg'


@pytest.mark.parametrize('error', [True, False])
def test_validation(dipole_map, error):
    assert utils.is_cea_map(dipole_map, error)
    assert utils.is_full_sun_synoptic_map(dipole_map, error)

import os
import numpy as np

import sunpy.map
import pytest

import pfsspy

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

def test_roll_to_CRLN(gong_map) :
    m = sunpy.map.Map(gong_map)

    # Test rolling raw map
    rolled = pfsspy.utils.roll_to_CRLN(m)
    assert rolled.data.shape == m.data.shape, "Rolled map shape has changed"
    assert rolled.meta["CRVAL1"] == 180., "CRLN0 not at left edge of map"

    # Test rolling map with non integer offset
    test_dat = m.data.copy()
    test_header = m.meta.copy()
    test_header["CRVAL1"] += 0.1
    rolled_non_int = pfsspy.utils.roll_to_CRLN(sunpy.map.Map(test_dat,test_header))
    assert rolled_non_int.data.shape == m.data.shape, "Rolled map shape has changed"
    assert rolled_non_int.meta["CRVAL1"] == 180., "CRLN0 not at left edge of map"

def test_interp_CAR2CEA(adapt_map) :
    m = pfsspy.utils.load_adapt(adapt_map)[0]
    interped = pfsspy.utils.interp_CAR2CEA(m)
    assert interped.meta['CTYPE1'] == "CRLN-CEA", "Not CEA Projection"
    assert interped.meta['CTYPE2'] == "CRLT-CEA", "Not CEA Projection"
    assert np.isclose(interped.meta['naxis2']*interped.meta['cdelt2']*np.pi/180,2.0), 'SinLat range is not 2.0'

def test_load_adapt(adapt_map) :
    adaptMapSequence = pfsspy.utils.load_adapt(adapt_map)
    assert isinstance(adaptMapSequence,sunpy.map.MapSequence)
    for map_ in adaptMapSequence : assert map_.meta['model'] == "ADAPT" 

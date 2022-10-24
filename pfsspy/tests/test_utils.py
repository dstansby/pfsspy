import astropy.units as u
import datetime
import numpy as np
import pytest
import sunpy.map

import pfsspy
from pfsspy import utils

from .example_maps import adapt_map, dipole_map, gong_map  # NoQA


def test_load_adapt(adapt_map):
    adaptMapSequence = utils.load_adapt(adapt_map)
    assert isinstance(adaptMapSequence, sunpy.map.MapSequence)
    for map_ in adaptMapSequence:
        assert map_.meta['model'] == "ADAPT"


def test_header_generation():
    ntheta = 180
    nphi = 90
    dtime = '2001-01-01 00:00:00'
    shape = [nphi, ntheta]
    header = pfsspy.utils.carr_cea_wcs_header(dtime, shape)
    assert header['LONPOLE'] == 0
    assert header['CTYPE1'] == 'CRLN-CEA'
    assert header['CTYPE2'] == 'CRLT-CEA'
    assert header['CDELT1'] == 360 / nphi
    np.testing.assert_almost_equal(
        header['CDELT2'], (180 / np.pi) * (2 / ntheta))

    assert header['CRPIX1'] == (nphi / 2) + 0.5
    assert header['CRPIX2'] == (ntheta / 2) + 0.5
    assert header['CRVAL1'] == 0
    assert header['CRVAL2'] == 0
    assert header['CUNIT1'] == 'deg'
    assert header['CUNIT2'] == 'deg'

    # Check corner coordinates are as expected
    data = np.random.rand(*shape)
    m = sunpy.map.Map(data.T, header)

    tols = {'rtol': 0, 'atol': 0.01 * u.deg}
    # Bottom left corner
    corner_coord = m.pixel_to_world(-0.5 * u.pix, -0.5 * u.pix)
    assert u.allclose(corner_coord.lat, -90 * u.deg, **tols)
    assert u.allclose(corner_coord.lon, 180 * u.deg, **tols)

    # Top right corner
    top_coord = m.pixel_to_world(89.5 * u.pix, 179.5 * u.pix)
    assert u.allclose(top_coord.lat, 90 * u.deg, **tols)
    assert u.allclose(corner_coord.lon, 180 * u.deg, **tols)


@pytest.mark.parametrize('error', [True, False])
def test_validation(dipole_map, error):
    assert utils.is_cea_map(dipole_map, error)
    assert utils.is_full_sun_synoptic_map(dipole_map, error)


def test_validation_not_full_map(dipole_map):
    dipole_map.meta['cdelt1'] = 0.001
    assert not utils.is_full_sun_synoptic_map(dipole_map)
    with pytest.raises(ValueError):
        utils.is_full_sun_synoptic_map(dipole_map, error=True)


def test_car_reproject(adapt_map):
    adapt_map = utils.load_adapt(adapt_map)[0]
    adapt_reproj = utils.car_to_cea(adapt_map)

    assert np.all(np.isfinite(adapt_map.data))
    assert np.all(np.isfinite(adapt_reproj.data))

    assert adapt_reproj.data.shape == adapt_map.data.shape
    for i in [1, 2]:
        assert adapt_reproj.meta[f'CTYPE{i}'][5:8] == 'CEA'

    with pytest.raises(ValueError, match='method must be one of'):
        utils.car_to_cea(adapt_map, method='gibberish')


def test_roll_map(gong_map):
    lh_edge_test = 0.0 * u.deg
    gong_map = sunpy.map.Map(gong_map)
    rolled_map = utils.roll_map(gong_map,
                                lh_edge_lon=lh_edge_test)

    # Test ref pixel rolled correctly
    # (-0.5, -0.5) is the bottom-left corner of the bottom-left pixel
    assert rolled_map.pixel_to_world(-0.5 * u.pixel,
                                     -0.5 * u.pixel).lon == lh_edge_test

    # Test output map is all finite
    assert np.all(np.isfinite(rolled_map.data))

    # Test output map is full sun synoptic
    assert utils.is_full_sun_synoptic_map(rolled_map, error=True)

    # Test reproject method error handling
    with pytest.raises(ValueError, match='method must be one of'):
        utils.roll_map(gong_map, method='gibberish')

    # Test left hand edge input type validation
    # 1. No Units
    with pytest.raises(TypeError,
                       match="has no 'unit' attribute"):
        utils.roll_map(adapt_map, lh_edge_lon=0)
    # 2. Incompatible units
    with pytest.raises(u.UnitsError,
                       match="must be in units convertible to 'deg'"):
        utils.roll_map(adapt_map, lh_edge_lon=0 * u.m)

    # Test left hand edge input range validation
    with pytest.raises(ValueError,
                       match='lh_edge_lon must be in'):
        utils.roll_map(adapt_map, lh_edge_lon=361 * u.deg)

def test_cea_header():
    # Assert default reference pixel is at 180 deg lon
    cea_default = utils.carr_cea_wcs_header(
        datetime.datetime(2020,1,1),
        [360,180]
    )
    assert cea_default['crval1'] == 0.0

    # Assert custom reference pixel is expected lon
    cea_shift = utils.carr_cea_wcs_header(
        datetime.datetime(2020,1,1),
        [360,180],
        map_center_longitude=10.0*u.deg
    )
    assert cea_shift['crval1'] == 10.0

    # Test reference pixel shift error handling
    # 1: No units
    with pytest.raises(u.UnitTypeError):
        cea_default = utils.carr_cea_wcs_header(
            datetime.datetime(2020,1,1),
            [360,180],
            map_center_longitude=0.0
        )
    # 2: Wrong Units
    with pytest.raises(u.UnitTypeError):
        cea_default = utils.carr_cea_wcs_header(
            datetime.datetime(2020,1,1),
            [360,180],
            map_center_longitude=0.0*u.m
        )
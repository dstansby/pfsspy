import astropy.constants as const
import astropy.coordinates as coord
import astropy.units as u
from astropy.tests.helper import quantity_allclose

import pytest

import matplotlib
import numpy as np
import pfsspy
import sunpy.map
import sunpy.util.exceptions

import pfsspy.coords
from pfsspy import tracing

from .example_maps import dipole_map, zero_map, dipole_result
matplotlib.use('Agg')

R_sun = const.R_sun


def test_expansion_factor(dipole_result):
    inp, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer()
    seed = coord.SkyCoord(0 * u.deg, 80 * u.deg, 1.1 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    assert field_line.expansion_factor > 1

    seed = coord.SkyCoord(0 * u.deg, -80 * u.deg, 1.1 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    assert field_line.expansion_factor > 1

    # This is a closed field line
    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.1 * R_sun, frame=out_frame)
    eq_field_line = tracer.trace(seed, out)[0]
    assert np.isnan(eq_field_line.expansion_factor)

    # Check that a field line near the equator has a bigger expansion
    # factor than one near the pole
    seed = coord.SkyCoord(0 * u.deg, 10 * u.deg, 2.49 * R_sun, frame=out_frame)
    pil_field_line = tracer.trace(seed, out)[0]
    assert pil_field_line.expansion_factor > field_line.expansion_factor


def test_field_line_polarity(dipole_result):
    input, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer()
    seed = coord.SkyCoord(0 * u.deg, 90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)
    assert field_line[0].polarity == 1

    seed = coord.SkyCoord(0 * u.deg, -90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)
    assert field_line[0].polarity == -1

    # This is a closed field line
    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.01 * R_sun, frame=out_frame)
    eq_field_line = tracer.trace(seed, out)[0]
    assert eq_field_line.polarity == 0


def test_footpoints(dipole_result):
    input, out = dipole_result
    out_frame = out.coordinate_frame

    tracer = tracing.PythonTracer(atol=1e-8, rtol=1e-8)

    def check_radius(coord, r):
        coord.representation_type = 'spherical'
        assert quantity_allclose(coord.radius, r)

    def check_open_fline(fline):
        check_radius(fline.solar_footpoint, const.R_sun)
        check_radius(fline.source_surface_footpoint, 2.5 * const.R_sun)

    seed = coord.SkyCoord(0 * u.deg, 90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_open_fline(field_line)

    seed = coord.SkyCoord(0 * u.deg, -90*u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_open_fline(field_line)

    seed = coord.SkyCoord(0 * u.deg, 0 * u.deg, 1.01 * R_sun, frame=out_frame)
    field_line = tracer.trace(seed, out)[0]
    check_radius(field_line.solar_footpoint, const.R_sun)
    check_radius(field_line.source_surface_footpoint, const.R_sun)


def test_shape(zero_map):
    # Test output map shapes
    input, out = zero_map
    nr = input.grid.nr
    nphi = input.grid.nphi
    ns = input.grid.ns

    out = pfsspy.pfss(input)
    alr, als, alp = out.al
    for comp in (alr, als, alp):
        assert np.all(comp == 0)

    assert alr.shape == (nphi + 1, ns + 1, nr)
    assert als.shape == (nphi + 1, ns, nr + 1)
    assert alp.shape == (nphi, ns + 1, nr + 1)

    br, bs, bp = out.bc
    for comp in (br, bs, bp):
        assert np.all(comp == 0)

    assert br.shape == (nphi, ns, nr + 1)
    assert bs.shape == (nphi, ns + 1, nr)
    assert bp.shape == (nphi + 1, ns, nr)

    bg = out.bg
    assert np.all(bg == 0)
    assert bg.shape == (nphi + 1, ns + 1, nr + 1, 3)


def test_input_output(dipole_result):
    # Smoke test of saving/loading files
    _, out = dipole_result
    out.save('test.npz')
    new_out = pfsspy.load_output('test.npz')
    assert (new_out.al[0] == out.al[0]).all()


def test_header_generation():
    dtime = None
    ntheta = 180
    nphi = 360
    shape = [nphi, ntheta]
    header = pfsspy.carr_cea_wcs_header(dtime, shape)
    assert header['LONPOLE'] == 0
    assert header['CTYPE1'] == 'CRLN-CEA'
    assert header['CTYPE2'] == 'CRLT-CEA'
    assert header['PV1_1'] == 1
    assert header['PV2_1'] == 1
    assert header['CDELT1'] == 360 / nphi
    np.testing.assert_almost_equal(
        header['CDELT2'], (180 / np.pi) * (2 / ntheta))

    assert header['CRPIX1'] == (nphi / 2) + 0.5
    assert header['CRPIX2'] == (ntheta / 2) + 0.5
    assert header['CRVAL1'] == 0
    assert header['CRVAL2'] == 0
    assert header['CUNIT1'] == 'deg'
    assert header['CUNIT2'] == 'deg'


def test_wrong_projection_error(dipole_map):
    dipole_map.meta['ctype1'] = 'HGLN-CAR'
    with pytest.raises(ValueError, match='must be CEA'):
        pfsspy.Input(dipole_map, 5, 2.5)


def test_non_map_input():
    with pytest.raises(ValueError, match='br must be a SunPy Map'):
        pfsspy.Input(np.random.rand(2, 2), 1, 1)

import numpy as np
import pfsspy
import pytest


@pytest.fixture
def zero_map():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))

    input = pfsspy.Input(br, nr, ns, nphi, rss)
    output = pfsspy.pfss(input)
    return input, output


def test_shape(zero_map):
    # Test output map shapes
    input, out = zero_map
    nr = input.nr
    nphi = input.np
    ns = input.ns

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

    assert br.shape == (nphi + 2, ns + 2, nr + 1)
    assert bs.shape == (nphi + 2, ns + 1, nr + 2)
    assert bp.shape == (nphi + 1, ns + 2, nr + 2)

    br, bs, bp = out.bg
    for comp in (br, bs, bp):
        assert np.all(comp == 0)

    assert br.shape == (nphi + 1, ns + 1, nr + 1)
    assert bs.shape == (nphi + 1, ns + 1, nr + 1)
    assert bp.shape == (nphi + 1, ns + 1, nr + 1)


def test_write_al(zero_map):
    _, output = zero_map
    # Test writing of vector potential
    output.save_a('a.netcdf')

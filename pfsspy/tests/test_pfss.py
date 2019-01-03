import numpy as np
import pfsspy


def test_zeros():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))

    input = pfsspy.Input(br, nr, ns, nphi, rss)

    out = pfsspy.pfss(input, output='a')
    for comp in (out.alr, out.als, out.alp):
        assert np.all(comp == 0)

    assert out.alr.shape == (nphi + 1, ns + 1, nr)
    assert out.als.shape == (nphi + 1, ns, nr + 1)
    assert out.alp.shape == (nphi, ns + 1, nr + 1)

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

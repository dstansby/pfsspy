import numpy as np
import pfsspy


def test_zeros():
    # Test a completely zero input
    ns = 30
    nphi = 20
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))
    output = 'none'

    out = pfsspy.pfss(br, nr, ns, nphi, rss, output='a')
    for comp in (out.alr, out.als, out.alp):
        assert np.all(comp == 0)

    assert out.alr.shape == (nphi + 1, ns + 1, nr)
    assert out.als.shape == (nphi + 1, ns, nr + 1)
    assert out.alp.shape == (nphi, ns + 1, nr + 1)

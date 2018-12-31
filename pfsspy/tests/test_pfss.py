import numpy as np
import pfsspy


def test_zeros():
    # Test a completely zero input
    ns = 30
    nphi = 30
    nr = 10
    rss = 2.5
    br = np.zeros((ns, nphi))
    output = 'none'

    out = pfsspy.pfss(br, nr, ns, nphi, rss, output='none')
    for comp in (out.alr, out.als, out.alp):
        assert np.all(comp == 0)

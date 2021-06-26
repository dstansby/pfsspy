import pytest 
from cartesian.utils import calc_cosine_heliocentric_angle, calc_mu, ExtentObj
import numpy as np 

def test_calc_mu():
    assert calc_mu(0,0)==1.0 

def test_calc_cosine_heliocentric_angle():
    coord1 = ExtentObj([[0,0],[200,200]],[200,200])
    coord2 = ExtentObj([[-200,-200],[0,0]],[200,200])
    mu1 = calc_cosine_heliocentric_angle(coord1)
    mu2 = np.flip(calc_cosine_heliocentric_angle(coord2))

    #Check 1: mu calculation should depend only on the radial coordinate, not on the sign.
    assert np.mean(np.abs(mu1-mu2))<1e-10

    # Check 2: Any input other than sunpy.Map or ExtentObj shouldn't work. 
    with pytest.raises(TypeError):
        calc_cosine_heliocentric_angle(np.array([0.0,0.0]))
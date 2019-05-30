import pytest
import numpy as np
from pypospack.eamtools import create_r

@pytest.mark.parametrize('rmax,rN',[(10.,1000)])
def test__create_r(rmax,rN):
    r = create_r(rmax,rN)
    assert isinstance(r,np.ndarray)
    assert r.max() == rmax
    assert r.size == rN




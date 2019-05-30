import pytest
import numpy as np
from pypospack.eamtools import create_rho

@pytest.mark.parametrize('rhomax,rhoN',[(10.,1000)])
def test__create_rho(rhomax,rhoN):
    rho = create_rho(rhomax,rhoN)
    assert isinstance(rho,np.ndarray)
    assert rho.max() == rhomax
    assert rho.size == rhoN




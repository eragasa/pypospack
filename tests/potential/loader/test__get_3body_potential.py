import pytest
from pypospack.potential import (
        ThreeBodyPotential,
        threebody_potential_names)
from pypospack.potential.loader import get_3body_potential

@pytest.mark.parametrize('name,symbols',[(k,['Ni']) for k in threebody_potential_names])
def test__get_pair_potential(name,symbols):
    o = get_3body_potential(name,symbols)
    assert isinstance(o, ThreeBodyPotential)

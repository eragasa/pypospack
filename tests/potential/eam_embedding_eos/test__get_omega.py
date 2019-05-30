import pytest

from pypospack.potentiam.eam_embedding_eos

@pytest.mark.parametrize(
        'a0','lattice_type',
        [
            (3.52,'fcc')
        ]
def test__get_omege(a0,lattice_type):
    omega = get_omega(a,lattice_type)


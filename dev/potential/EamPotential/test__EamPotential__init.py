import pytest
from collections import OrderedDict
from pypospack.potential import (PairPotential,
                                 EamDensityFunction,
                                 EamEmbeddingFunction,
                                 EamEmbeddingEquationOfState,
                                 pair_potential_names,
                                 eam_density_names,
                                 eam_embedding_names
                                 )
from eam import EamPotential

@pytest.mark.parametrize('func_pair_name',[(k) for k in pair_potential_names])
@pytest.mark.parametrize('func_density_name',[(k) for k in eam_density_names])
@pytest.mark.parametrize('func_embedding_name',[(k) for k in eam_embedding_names])
def test____init__1sym(func_pair_name,func_density_name,func_embedding_name):
    symbols = ['Ni']

    eam = EamPotential(
            symbols=symbols,
            func_pair=func_pair_name,
            func_density=func_density_name,
            func_embedding=func_embedding_name)
    assert isinstance(eam.obj_pair,PairPotential)
    assert isinstance(eam.obj_density,EamDensityFunction)
    assert isinstance(eam.obj_embedding,EamEmbeddingFunction)

    #<------ testing eam.symbols
    assert type(eam.symbols) is list
    assert eam.symbols == symbols


import pytest
import inspect

from pypospack.potential.pair_general_lj import func_pair_generalized_lj_w_cutoff
from pypospack.potential import GeneralizedLennardJonesPotential

def test__GeneralizedLennardJonesPotential__potential_type():
    potential_type = 'general_lj'
    assert GeneralizedLennardJonesPotential.potential_type == potential_type

def test__GeneralizedLennardJonesPotential__pair_potential_parameters():
    pair_potential_parameters = ['b1','b2','r1','V0','delta','rc','hc','h0']
    assert GeneralizedLennardJonesPotential.pair_potential_parameters \
            == pair_potential_parameters

def test__GeneralizedLennardJonesPotential__callable_func():
    assert GeneralizedLennardJonesPotential.callable_func \
            == func_pair_generalized_lj_w_cutoff

def test____init__():
    symbols = ['Ni']
    potential_type = 'general_lj'
    pair_potential_parameters = ['b1','b2','r1','V0','delta','rc','hc','h0']

    o = GeneralizedLennardJonesPotential(symbols=symbols)
    assert o.potential_type == GeneralizedLennardJonesPotential.potential_type
    assert o.pair_potential_parameters == pair_potential_parameters
  

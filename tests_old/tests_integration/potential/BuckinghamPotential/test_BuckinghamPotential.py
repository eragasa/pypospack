import pytest
from collections import OrderedDict
import numpy as np
import pypospack.potential as potential

def test__import__pypospack_potential():
    from pypospack.potential import BuckinghamPotential

def test__import__pypospack_potentials_morse():
    from pypospack.potential import BuckinghamPotential

def test__2element____init__():
    symbols = ['Mg','O']

    try:
        buck = potential.BuckinghamPotential(symbols=symbols)
    except:
        pytest.fail()

    #<--- testing variables
    potential_type = 'buckingham'
    symbol_pairs = [['Mg','Mg'],['Mg','O'],['O','O']]
    parameter_names = ['chrg_Mg','chrg_O',
            'MgMg_A','MgMg_rho','MgMg_C',
            'MgO_A','MgO_rho','MgO_C',
            'OO_A','OO_rho','OO_C']

    #<--- test attribute potential_type
    assert type(buck.potential_type) is str
    assert buck.potential_type == 'buckingham'

    #<--- test attribute symbol_pairs
    assert type(buck.symbol_pairs) is list
    assert len(buck.symbol_pairs) == len(symbol_pairs)
    assert buck.symbol_pairs == symbol_pairs

    #<--- test attribute parameter_names
    assert type(buck.parameter_names) is list
    assert len(buck.parameter_names) == len(parameter_names)
    assert buck.parameter_names == parameter_names

    #<--- test attribute parameters
    assert type(buck.parameters) is OrderedDict
    for name,value in buck.parameters.items():
        assert value is None
    for name in buck.parameter_names:
        assert name in buck.parameter_names

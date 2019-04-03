"""
This module tests the pypospack.potential.MorsePotential
"""
from collections import OrderedDict
import pytest
import numpy as np

from pypospack.potential import PairPotential
from pypospack.potential import GeneralizedLennardJonesPotential

def get_testing_sets(testing_name=None):
    testing_sets = OrderedDict()
    testing_sets['Ni'] = OrderedDict()
    testing_sets['Ni']['symbols'] = ['Ni']
    testing_sets['Ni']['symbol_pairs'] = [['Ni','Ni']]
    testing_sets['Ni']['parameter_names'] = [
        'NiNi_b1', 'NiNi_b2', 'NiNi_r1', 'NiNi_V0', 'NiNi_delta'
    ]
    if testing_name is None:
        return [k for k in testing_sets], [k for k in testing_sets.values()]
    else:
        return testing_sets[testing_name]
testing_names, testing_sets = get_testing_sets()

def dev____init__():
    testing_set = get_testing_sets(testing_name='Ni')

    o = GeneralizedLennardJonesPotential(symbols=testing_set['symbols'])

    print('o.symbol_pairs:{}'.format(o.symbol_pairs))
    print('o.parameter_names:{}'.format(o.parameter_names))
    for k,v in o.parameters.items():
        print("{}:{}".format(k,v))


@pytest.mark.parametrize("testing_set",testing_sets,ids=testing_names)
def test____init__(testing_set):
    o = GeneralizedLennardJonesPotential(symbols=testing_set['symbols'])

    assert isinstance(o,PairPotential)
    assert isinstance(o.symbols,list)
    assert isinstance(o.symbol_pairs,list)
    assert isinstance(o.parameter_names,list)
    assert isinstance(o.parameters,OrderedDict)

    assert o.symbols == testing_set['symbols']
    assert o.symbol_pairs == testing_set['symbol_pairs']
    assert o.parameter_names == testing_set['parameter_names']
    assert set(o.parameters.keys()) == set(testing_set['parameters'].keys())
    for k,v in testing_set['parameters'].items():
        o.parameters[k] = v


if __name__ == "__main__":
    dev____init__()

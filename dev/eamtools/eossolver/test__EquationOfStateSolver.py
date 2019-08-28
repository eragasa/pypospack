import pytest
from collections import OrderedDict

import numpy as np

from pypospack.eamtools import create_r,create_rho
from eossolver import EquationOfStateSolver
args = ['symbols']
testcase_dict = OrderedDict()
testcase_dict['Ni'] = OrderedDict()
testcase_dict['Ni']['symbols'] = ['Ni']

testcase_names =  [k for k in testcase_dict.keys()]
testcases = [v for v in testcase_dict.values()]
testcase_values = [tuple(v.values()) for v in testcases]

from eossolver import determine_pair_names
@pytest.mark.parametrize(
    "symbols,exp_pair_names",
    [
        (
            ['Ni'],
            ['NiNi']
        ),(
            ['Ni','Al'],
            ['NiNi','AlNi','AlAl']
        )
    ]
)
def test__determine_pair_names(symbols,exp_pair_names):
    pair_names = determine_pair_names(symbols)
    assert pair_names == exp_pair_names

@pytest.mark.parametrize(
    "symbols,pair_names",
    [
        (
            ['Ni'],
            ['NiNi'])
    ]
)
def test__EquationOfStateSolver(symbols,pair_names):
    rmax = 10
    rN = 1000
    r = create_r(rmax,rN)

    rhomax = 10.
    rhoN = 1000
    rho = create_rho(rhomax,rN)

    o = EquationOfStateSolver(symbols,r,rho)
    assert o.symbols == symbols
    assert np.array_equal(o.r,r)
    assert np.array_equal(o.rho,rho)
    
    assert isinstance(o.potentials,OrderedDict)
    for k in EquationOfStateSolver.potential_types:
        assert k in o.potentials
    for pn in pair_names:
        assert pn in o.potentials['pair']
        assert isinstance(o.potentials['pair'][pn],OrderedDict)
        assert 'formalism' in o.potentials['pair'][pn]
        assert 'parameters' in o.potentials['pair'][pn]
        assert 'evaluation' in o.potentials['pair'][pn]
    for s in symbols:
        assert s in o.potentials['density']
        assert isinstance(o.potentials['density'][s],OrderedDict)
        assert 'formalism' in o.potentials['density'][s]
        assert 'parameters' in o.potentials['density'][s]
        assert 'evaluation' in o.potentials['density'][s]
    for s in symbols:
        assert s in o.potentials['embedding']
        assert isinstance(o.potentials['embedding'][s],OrderedDict)
        assert 'formalism' in o.potentials['embedding'][s]
        assert 'parameters' in o.potentials['embedding'][s]
        assert 'evaluation' in o.potentials['embedding'][s]

import pytest
import os
from collections import OrderedDict

import pypospack.utils
from pypospack.eamtools import SeatonSetflReader
from potentialfitter import determine_eam_setfl_pairs
from potentialfitter import EamPotentialFitter

## testing functions
from pypospack.potential.eamdens_mishin2003 import func_mishin2003_density_w_cutoff
from pypospack.potential.pair_general_lj import func_generalized_lj_w_cutoff

# Mishin2003
potentials = OrderedDict()
potentials['setfl_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
potentials['density'] = OrderedDict()
potentials['density']['Ni'] = OrderedDict()
potentials['density']['Ni']['formalism'] = func_mishin2003_density_w_cutoff
potentials['density']['Ni']['param'] = OrderedDict()
potentials['density']['Ni']['param']['r0'] = -3.138
potentials['density']['Ni']['param']['A0'] = 1.
potentials['density']['Ni']['param']['B0'] = 1.1914e4
potentials['density']['Ni']['param']['C0'] = 2.0329e2
potentials['density']['Ni']['param']['y'] = 1.9521
potentials['density']['Ni']['param']['gamma'] = 1.6802e3
potentials['density']['Ni']['param']['rc']=5.168
potentials['density']['Ni']['param']['h']=3.32

potentials['pair'] = OrderedDict()
potentials['pair']['NiNi'] = OrderedDict()
potentials['pair']['NiNi']['formalism'] = func_generalized_lj_w_cutoff 
potentials['pair']['NiNi']['param'] = OrderedDict()
potentials['pair']['NiNi']['param']['b1'] = 4.7067e-3     # no units
potentials['pair']['NiNi']['param']['b2'] = 0.15106       # no units
potentials['pair']['NiNi']['param']['r1'] = 3.8673e-4      # angs
potentials['pair']['NiNi']['param']['delta'] = 3.6046e3   # eV
potentials['pair']['NiNi']['param']['V0'] = -3.5126e3     # eV
potentials['pair']['NiNi']['param']['rc'] = 5.168         # angs
potentials['pair']['NiNi']['param']['h'] = 3.3228         # angs

@pytest.mark.parametrize('symbols',[
    (['Ni']),
    (['Ni','Al'])])
def test__EamPotentialFitter__init__(symbols):
    o = EamPotentialFitter(symbols)
    assert isinstance(o.pair_potentials,OrderedDict)
    assert isinstance(o.density_functions,OrderedDict)
    assert isinstance(o.embedding_functions,OrderedDict)

    assert list(o.pair_potentials.keys()) \
            == list(["".join(k) for k in determine_eam_setfl_pairs(symbols)])
    assert list(o.density_functions.keys()) == symbols
    assert list(o.embedding_functions.keys()) == symbols

    for k1 in ['p0','popt']:
        assert k1 in o.parameters
        
        for k2 in ['pair','embedding','density']:
            assert k2 in o.parameters[k1]

        for pair in determine_eam_setfl_pairs(symbols):
            assert "".join(pair) in o.parameters[k1]['pair']
            assert o.parameters[k1]['pair']["".join(pair)] is None

        for s in symbols:
            assert s in o.parameters[k1]['embedding']
            assert o.parameters[k1]['embedding'][s] is None

        for s in symbols:
            assert s in o.parameters[k1]['density']
            assert o.parameters[k1]['embedding'][s] is None

    for k in ['pair','embedding','density']:
        assert isinstance(o.formalisms[k],OrderedDict)

    for pair in determine_eam_setfl_pairs(symbols):
        pair_name  = "".join(pair)
        assert o.formalisms['pair'][pair_name] is None
    
    for s in symbols:
        assert o.formalisms['density'][s] is None

    for s in symbols:
        assert o.formalisms['embedding'][s] is None


    assert all([v is None for v in o.pair_potentials.values()])
    assert all([v is None for v in o.density_functions.values()])
    assert all([v is None for v in o.embedding_functions.values()])

    assert o.setfl_reader is None

@pytest.mark.parametrize('symbols,setfl_filename',[
        (
            ['Ni','Al'],
            os.path.join(
                pypospack.utils.get_pypospack_root_directory(),
                'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
        )
    ])
def test__read_setfl_filename(symbols,setfl_filename):
    o = EamPotentialFitter(symbols)
    o.read_setfl_file(filename=setfl_filename)

    assert isinstance(o.setfl_reader,SeatonSetflReader)

@pytest.mark.parametrize(
        'symbols,setfl_fn,func_pair_potential,symbol_pair,param0',
        [
            (  
                ['Ni','Al'],
                potentials['setfl_fn'],
                potentials['pair']['NiNi']['formalism'],
                ['Ni','Ni'],
                potentials['pair']['NiNi']['param'],
            )
        ]
)
def test__fit_potential_pair(symbols,setfl_fn,func_pair_potential,symbol_pair,param0):
    o = EamPotentialFitter(symbols)
    o.read_setfl_file(filename=setfl_fn)

    o.fit_potential_pair(
            func_pair_potential=func_pair_potential,
            symbol_pair=symbol_pair,
            param0=param0,
            rlow=1.5)
if __name__ == "__main__":
    setfl_fn = os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
    symbols = ['Ni']
    o = EamPotentialFitter(symbols)
    o.read_setfl_file(filename=setfl_fn)
    o.fit_potential_pair(
            func_pair_potential=potentials['pair']['NiNi']['formalism'],
            symbol_pair=['Ni','Ni'],
            param0=potentials['pair']['NiNi']['param'],
            rlow=2.)

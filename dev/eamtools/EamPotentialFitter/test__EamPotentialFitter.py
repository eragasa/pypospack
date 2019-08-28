import pytest
import os
from collections import OrderedDict
import pypospack.utils
from pypospack.potential.pair_general_lj import func_pair_generalized_lj_w_cutoff
from potentialfitter import EamPotentialFitter

test_cases = OrderedDict()
test_cases['punmishin2009'] = OrderedDict()
test_cases['punmishin2009']['symbols'] = ['Ni']
test_cases['punmishin2009']['setfl_fn'] = os.path.join(
        pypospack.utils.get_pypospack_root_directory(),
        'data','potentials','Ni__eam','Mishin-Ni-Al-2009.eam.alloy')
test_cases['punmishin2009']['pair'] = OrderedDict()
test_cases['punmishin2009']['pair']['NiNi'] = OrderedDict()
test_cases['punmishin2009']['pair']['NiNi']['formalism'] = func_pair_generalized_lj_w_cutoff
test_cases['punmishin2009']['pair']['NiNi']['pair'] = ['Ni', 'Ni']
test_cases['punmishin2009']['pair']['NiNi']['param'] = OrderedDict()
test_cases['punmishin2009']['pair']['NiNi']['param']['b1'] = 4.7067e-3     # no units
test_cases['punmishin2009']['pair']['NiNi']['param']['b2'] = 0.15106       # no units
test_cases['punmishin2009']['pair']['NiNi']['param']['r1'] = 3.8673e-4      # angs
test_cases['punmishin2009']['pair']['NiNi']['param']['delta'] = 3.6046e3   # eV
test_cases['punmishin2009']['pair']['NiNi']['param']['V0'] = -3.5126e3     # eV
test_cases['punmishin2009']['pair']['NiNi']['param']['rc'] = 5.168         # angs
test_cases['punmishin2009']['pair']['NiNi']['param']['hc'] = 3.3228         # angs
test_cases['punmishin2009']['pair']['NiNi']['param']['h0'] = 3.3228         # angs


@pytest.mark.parametrize('symbols',
                         [ (test_cases['punmishin2009']['symbols']) ])
def test____init__(symbols):
    o = EamPotentialFitter(symbols)

@pytest.mark.parametrize('symbols,setfl_fn',
                         [ (test_cases['punmishin2009']['symbols'],
                            test_cases['punmishin2009']['setfl_fn']) ])
def test__read_setfl_file(symbols, setfl_fn):
    o = EamPotentialFitter(symbols)
    o.read_setfl_file(filename=setfl_fn)


@pytest.mark.parametrize('symbols,setfl_fn,pair_definition',
                         [ (test_cases['punmishin2009']['symbols'],
                            test_cases['punmishin2009']['setfl_fn'],
                            test_cases['punmishin2009']['pair']['NiNi']) ])
def test__fit_potential_pair(symbols, setfl_fn, pair_definition):

    o = EamPotentialFitter(symbols)
    o.read_setfl_file(filename=setfl_fn)

    o.fit_potential_pair(
            func_pair_potential=pair_definition['formalism'],
            symbol_pair=pair_definition['pair'],
            param0=pair_definition['param'])

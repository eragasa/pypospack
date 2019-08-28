import pytest

import inspect
import numpy as np
from pypospack.eamtools import create_r
from pypospack.potential.pair_general_lj import (func_cutoff_mishin2004,
                                                 func_pair_generalized_lj,
                                                 func_pair_generalized_lj_w_cutoff)

parameters = {
        'b1':4.7067e-3,
        'b2':0.15106,
        'r1':3.8663e-4,
        'V0':-3.5126e3,
        'delta':3.6046e3,
        'rc':5.168,
        'hc':3.323,
        'h0':1.500
    }

r = create_r(6.0,1000)

def test__func_cutoff_mishin2004__r_numpy():
    param = [parameters[v] for v in inspect.signature(func_cutoff_mishin2004).parameters if v != 'r']
    psi = func_cutoff_mishin2004(r, *param)
    assert isinstance(psi, np.ndarray)

def test__func_pair_generalized_lj__r_numpy():
    param = [parameters[v] for v in inspect.signature(func_pair_generalized_lj).parameters if v != 'r']
    phi = func_pair_generalized_lj(r,*param)
    assert isinstance(phi, np.ndarray)

def test__test__func_pair_generalized_lj__r_numpy():
    phi = func_pair_generalized_lj_w_cutoff(r, **parameters)
    assert isinstance(phi, np.ndarray)


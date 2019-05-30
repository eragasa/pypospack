import pytest
import numpy as np
from collections import OrderedDict
from pypospack.potential import LennardJonesPotential

testcase = OrderedDict()
testcase['symbols'] = ['Ar']
testcase['parameters'] = OrderedDict()
testcase['parameters']['ArAr_sigma'] = 1.0
testcase['parameters']['ArAr_epsilon'] = 1.0
testcase['parameter_names'] = ['r_cut_global', 'ArAr_epsilon', 'ArAr_sigma', 'ArAr_r_cut_pair', 'ArAr_r_cut_coulomb']

def test__static_variable__global_parameters():
    assert isinstance(LennardJonesPotential.global_potential_parameters,list)
    assert LennardJonesPotential.global_potential_parameters == ['r_cut_global']

def test__static_variable__pair_parameters():
    assert isinstance(LennardJonesPotential.pair_potential_parameters,list)
    assert LennardJonesPotential.pair_potential_parameters == ['epsilon','sigma','r_cut_pair','r_cut_coulomb']

def test____init__1sym():
    o = LennardJonesPotential(symbols=testcase['symbols'])

    assert isinstance(o.parameter_names,list)
    assert o.parameter_names == testcase['parameter_names']
    
    assert isinstance(o.parameters,OrderedDict)
    assert len(o.parameters) == len(o.parameter_names)
    assert all([v in o.parameters.keys() for v in testcase['parameter_names']])

def test__evaluate():
    r = np.linspace(1,1000,1000)/1000*3
    o = LennardJonesPotential(symbols=testcase['symbols'])
    o.evaluate(
            r=r,
            parameters=testcase['parameters'])
    assert isinstance(o.parameters,OrderedDict)
    assert isinstance(o.potential_evaluations,OrderedDict)
    assert all([isinstance(v,np.ndarray) for k,v in o.potential_evaluations.items()])
    assert all([v.shape == r.shape for k,v in o.potential_evaluations.items()])

def dev____init____1sym():
    o = LennardJonesPotential(symbols=testcase['symbols'])
    print(o.parameter_names)
    print(o.parameters.keys())

def dev__evaluate():
    
    r = np.linspace(1,1000,1000)/1000*3
    o = LennardJonesPotential(symbols=testcase['symbols'])
    o.evaluate(
            r=r,
            parameters=testcase['parameters'])
    print(o.parameters)
    print(o.potential_evaluations)

if __name__ == '__main__':
    dev____init____1sym()
    dev__evaluate()

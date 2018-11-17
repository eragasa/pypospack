import pytest
from collections import OrderedDict
from pypospack.potential import EamDensityFunction
from pypospack.potential import ExponentialDensityFunction

symbols = ['Al']
a0 = 4.046
latt_type = 'fcc'
expected_parameter_names = [
        '{}_rho0'.format(symbols[0]),
        '{}_beta'.format(symbols[0]),
        '{}_r0'.format(symbols[0])
        ]

import numpy as np
rcut = 10
rho = np.linspace(0,rcut,100)

def get_expected_parameter_names(symbols):
    parameter_names = []
    for s1 in symbols:
        pass
        

def test____init():
    dens = ExponentialDensityFunction(symbols=symbols)
    
    assert isinstance(dens,EamDensityFunction)

    assert type(dens.symbols) is list
    assert len(dens.symbols) == len(symbols)
    assert dens.symbols == symbols

    assert type(dens.parameter_names) is list
    assert len(dens.parameter_names) == len(expected_parameter_names)
    assert dens.parameter_names == expected_parameter_names

    assert type(dens.parameters) is OrderedDict
    assert len(dens.parameters) is len(dens.parameter_names)
    for p in dens.parameter_names:
        assert p in dens.parameters

def test__evaluate():
    dens = ExponentialDensityFunction(symbols=symbols)

    test_parameters = OrderedDict([
        ('Al_rho0',1.0),
        ('Al_beta',1.0),
        ('Al_r0',1.0)])

    dens.evaluate(r=rho,parameters=test_parameters)

def test__evaluation__with_bad_parameter_name():
    dens = ExponentialDensityFunction(symbols=symbols)

    test_parameters = OrderedDict([
        ('Al_rho',1.0),
        ('Al_beta',1.0),
        ('Al_r0',1.0)])

    with pytest.raises(KeyError):
        dens.evaluate(r=rho,parameters=test_parameters)

if __name__ == "__main__":
    pass

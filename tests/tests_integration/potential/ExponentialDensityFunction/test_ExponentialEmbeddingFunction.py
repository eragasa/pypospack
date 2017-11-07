import pytest
from collections import OrderedDict
from pypospack.potential import EamDensityFunction

def test__import__from_pypospack_potential():
    from pypospack.potential import ExponentialDensityFunction

def test_1sym____init__():
    #<--- variables
    symbols = ['Ni']

    #<--- setup
    from pypospack.potential import ExponentialDensityFunction

    #<--- code being tested
    dens = ExponentialDensityFunction(symbols=symbols)

    #<--- expected results
    parameter_names = ['Ni_rho0','Ni_beta','Ni_r0']
    
    #<--- testing expected results
    assert isinstance(dens,EamDensityFunction)

    assert type(dens.symbols) is list
    assert len(dens.symbols) == len(symbols)
    assert dens.symbols == symbols

    assert type(dens.parameter_names) is list
    assert len(dens.parameter_names) == len(parameter_names)
    assert dens.parameter_names == parameter_names

    assert type(dens.parameters) is OrderedDict
    assert len(dens.parameters) is len(dens.parameter_names)
    for p in dens.parameter_names:
        assert p in dens.parameters

    assert dens.density is None

def test_2sym__init__():
    #<--- variables
    symbols = ['Ni','Al']

    #<--- setup
    from pypospack.potential import ExponentialDensityFunction

    #<--- code being tested
    dens = ExponentialDensityFunction(symbols=symbols)

    #<--- expected results
    parameter_names = ['Ni_rho0','Ni_beta','Ni_r0',
                       'Al_rho0','Al_beta','Al_r0']
    
    #<--- testing expected results
    assert isinstance(dens,EamDensityFunction)

    assert type(dens.symbols) is list
    assert len(dens.symbols) == len(symbols)
    assert dens.symbols == symbols

    assert type(dens.parameter_names) is list
    assert len(dens.parameter_names) == len(parameter_names)
    assert dens.parameter_names == parameter_names

    assert type(dens.parameters) is OrderedDict
    assert len(dens.parameters) is len(dens.parameter_names)
    for p in dens.parameter_names:
        assert p in dens.parameters

    assert dens.density is None

if __name__ == "__main__":
    from pypospack.potential import ExponentialDensityFunction
    symbols = ['Ni']
    dens = ExponentialDensityFunction(symbols=symbols)
    for p in dens.parameters:
        print(p)
    
    symbols = ['Ni','Al']
    dens = ExponentialDensityFunction(symbols=symbols)
    for p in dens.parameters:
        print(p)

    
    symbols = ['Ni','Al']
    dens = ExponentialDensityFunction(symbols=symbols)
    for p in dens.parameters:
        print(p)
    symbols = ['Ni','Al']
    dens = ExponentialDensityFunction(symbols=symbols)
    for p in dens.parameters:
        print(p)

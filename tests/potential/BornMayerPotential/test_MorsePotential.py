import pytest
from collections import OrderedDict
import numpy as np
import pypospack.potential as potential

def test__import__pypospack_potential():
    from pypospack.potential import MorsePotential

def test__import__pypospack_potentials_morse():
    from pypospack.potentials.morse import MorsePotential

def test__1element____init__():
    symbols = ['Ni']
    symbol_pairs = [['Ni','Ni']]
    parameter_names = ['NiNi_D0','NiNi_a','NiNi_r0']
    parameters = OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813

    r_max = 11.
    N_r = 500

    try:
        morse = potential.MorsePotential(symbols=symbols)
    except:
        pytest.fail()

    #<-- test attribute symbols
    assert type(morse.symbols) is list
    assert len(morse.symbols) == len(symbols)
    assert morse.symbols == symbols

    #<--- test attribute symbol_pairs
    assert type(morse.symbol_pairs) is list
    assert len(morse.symbol_pairs) == len(symbol_pairs)
    assert morse.symbol_pairs == symbol_pairs

    #<-- test attribute parameter_names
    assert type(morse.parameter_names) is list
    assert len(morse.parameter_names) == len(parameter_names)
    assert morse.parameter_names == parameter_names
    #<-- test attribute parameters
    assert type(morse.parameters) is OrderedDict
    for name,value in morse.parameters.items():
        assert value is None
    assert len(morse.parameters) == len(morse.parameter_names)
    for name in morse.parameter_names:
        assert name in morse.parameters

def test__1element__evaluate():
    symbols = ['Ni']
    symbol_pairs = [['Ni','Ni']]
    parameter_names = ['NiNi_D0','NiNi_a','NiNi_r0']
    parameters = OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813

    r_max = 11.
    N_r = 500
    r = r_max * np.linspace(1,100,N_r)/100

    try:
        morse = potential.MorsePotential(symbols=symbols)
        morse.evaluate(r,parameters)
    except:
        pytest.fail()

    assert isinstance(morse.potential, OrderedDict)
    for pair_key,pot in morse.potential.items():
        assert isinstance(pot,np.ndarray)
        assert pot.shape == r.shape

def test__2element___init__():
    symbols = ['Ni','Al']
    parameters= OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813
    parameters['NiAl_D0'] = 0.001114
    parameters['NiAl_a'] = 3.429506
    parameters['NiAl_r0'] = 2.6813
    parameters['AlAl_D0'] = 0.001114
    parameters['AlAl_a'] = 3.429506
    parameters['AlAl_r0'] = 2.6813

    symbol_pairs = [['Ni','Ni'],['Ni','Al'],['Al','Al']]
    parameter_names = ['NiNi_D0','NiNi_a','NiNi_r0',
                       'NiAl_D0','NiAl_a','NiAl_r0',
                       'AlAl_D0','AlAl_a','AlAl_r0']

    try:
        morse = potential.MorsePotential(symbols=symbols)
    except:
        pytest.fail()

    assert type(morse.symbols) is list
    assert type(morse.symbol_pairs) is list
    assert morse.symbol_pairs == symbol_pairs
    assert type(morse.parameter_names) is list
    assert morse.parameter_names == parameter_names
if __name__ == "__main__":

    symbols = ['Ni']
    symbol_pairs = [['Ni','Ni']]
    parameter_names = ['NiNi_D0','NiNi_a','NiNi_r0']
    parameters = OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813

    r_max = 11.
    N_r = 500
    r = r_max * np.linspace(1,100,N_r)/100

    try:
        morse = potential.MorsePotential(symbols=symbols)
        morse.evaluate(r,parameters)
    except:
        print(morse.parameters)

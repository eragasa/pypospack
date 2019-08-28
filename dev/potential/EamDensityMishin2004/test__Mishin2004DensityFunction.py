from pypospack.potential import EamDensityFunction
from eamdens_mishin2004 import Mishin2004DensityFunction

def test__Mishin2004DensityFunction__potential_type():
    assert Mishin2004DensityFunction.potential_type == 'eamdens_mishin2004'

def test__Mishin2004DensityFunction__density_function_parameters():
    expected_names = ['A0','B0','C0','y','gamma', 'rc', 'hr', 'h0']
    assert Mishin2004DensityFunction.density_function_parameters == expected_names

def test____init__():
    symbols = ['Ni']
    o = Mishin2004DensityFunction(symbols=symbols)
    assert isinstance(o, EamDensityFunction)
    assert o.potential_type == Mishin2004DensityFunction.potential_type
    assert isinstance(o.parameter_names, list)

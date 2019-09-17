from pypospack.potential import EamDensityFunction
from pypospack.potential import Mishin2004DensityFunction

def test__Mishin2004DensityFunction__potential_type():
    assert Mishin2004DensityFunction.potential_type == 'eamdens_mishin2004'

def test__Mishin2004DensityFunction__density_function_parameters():
    expected_names = ['r0','A0','B0','C0','y','gamma', 'rc', 'hc', 'h0']
    assert Mishin2004DensityFunction.density_function_parameters == expected_names

def test____init__():
    symbols = ['Ni']
    o = Mishin2004DensityFunction(symbols=symbols)
    assert isinstance(o, EamDensityFunction)
    assert o.potential_type == Mishin2004DensityFunction.potential_type
    assert isinstance(o.parameter_names, list)

def test__evaluate():
    parameters = {
            'Ni_r0':-4.74407425932,
            'Ni_A0':0.547690751136,
            'Ni_B0':9209.25239146,
            'Ni_C0':74.9843643,
            'Ni_gamma':1523.14810285,
            'Ni_y':2.23169424406,
            'Ni_rc':5.76416508323,
            'Ni_hc':0.305546154972,
            'Ni_h0':1.5262573968}
    r=3.5
    symbols = ['Ni']
    o = Mishin2004DensityFunction(symbols=symbols)
    o.evaluate(r,parameters)

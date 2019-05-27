"""
This module tests the pypospack.potential.MorsePotential
"""
from collections import OrderedDict
import pytest
import numpy as np

testing_set = OrderedDict()
testing_set['Ni'] = OrderedDict()
testing_set['Ni']['symbols'] = ['Ni']
testing_set['Ni']['symbol_pairs'] = [['Ni', 'Ni']]
testing_set['Ni']['parameter_names'] = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
testing_set['Ni']['parameters'] = OrderedDict()
testing_set['Ni']['parameters']['NiNi_D0'] = 0.001114
testing_set['Ni']['parameters']['NiNi_a'] = 3.429506
testing_set['Ni']['parameters']['NiNi_r0'] = 2.6813
testing_set['NiAl'] = OrderedDict()
testing_set['NiAl']['symbols'] = ['Ni', 'Al']
testing_set['NiAl']['symbol_pairs'] = [['Ni', 'Ni'], ['Ni', 'Al'], ['Al', 'Al']]
testing_set['NiAl']['parameters'] = OrderedDict()
testing_set['NiAl']['parameters']['NiNi_D0'] = 0.001114
testing_set['NiAl']['parameters']['NiNi_a'] = 3.429506
testing_set['NiAl']['parameters']['NiNi_r0'] = 2.6813
testing_set['NiAl']['parameters']['NiAl_D0'] = 0.001114
testing_set['NiAl']['parameters']['NiAl_a'] = 3.429506
testing_set['NiAl']['parameters']['NiAl_r0'] = 2.6813
testing_set['NiAl']['parameters']['AlAl_D0'] = 0.001114
testing_set['NiAl']['parameters']['AlAl_a'] = 3.429506
testing_set['NiAl']['parameters']['AlAl_r0'] = 2.6813


Ni_testing_set = OrderedDict()
Ni_testing_set['symbols'] = ['Ni']
Ni_testing_set['symbol_pairs'] = [['Ni', 'Ni']]
Ni_testing_set['parameter_names'] = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
Ni_testing_set['parameters'] = OrderedDict([
    ('NiNi_D0', 0.001114),
    ('NiNi_a', 3.429506),
    ('NiNi_r0', 2.6813)
])
Ni_testing_set['r_max'] = 11.0
Ni_testing_set['N_r'] = 500

NiAl_testing_set = OrderedDict()
NiAl_testing_set['symbols'] = ['Ni', 'Al']
NiAl_testing_set['symbol_pairs'] = [['Ni', 'Ni'], ['Ni', 'Al'], ['Al', 'Al']]
NiAl_testing_set['parameters'] = OrderedDict([
    ('NiNi_D0', 0.001114), ('NiNi_a', 3.429506), ('NiNi_r0', 2.6813),
    ('NiAl_D0', 0.001114), ('NiAl_a', 3.429506), ('NiAl_r0', 2.6813),
    ('AlAl_D0', 0.001114), ('AlAl_a', 3.429506), ('AlAl_r0', 2.6813)
])
NiAl_testing_set['parameter_names'] = [
    'NiNi_D0', 'NiNi_a', 'NiNi_r0',
    'NiAl_D0', 'NiAl_a', 'NiAl_r0',
    'AlAl_D0', 'AlAl_a', 'AlAl_r0'
]
NiAl_testing_set['r_max'] = 11.0
NiAl_testing_set['N_r'] = 500

testing_sets = [(Ni_testing_set), (NiAl_testing_set)]
testing_names = ['Ni', 'NiAl']

class TestMorsePotential(object):

    def test_initialize_no_args(self):
        from pypospack.potential import MorsePotential

        with pytest.raises(Exception) as e:
            p = MorsePotential()

    @pytest.mark.parametrize("testing_set", testing_sets, ids=testing_names)
    def test_initialize(self, testing_set):

        print('symbols:', testing_set['symbols'])
        from pypospack.potential import PairPotential
        from pypospack.potential import MorsePotential
        o = MorsePotential(symbols=testing_set['symbols'])
        
        assert isinstance(o,PairPotential)
        assert isinstance(o.symbols,list)
        assert isinstance(o.symbol_pairs,list)
        assert isinstance(o.parameter_names,list)
        assert isinstance(o.parameters,OrderedDict)

        assert o.symbols == testing_set['symbols']
        assert o.symbol_pairs == testing_set['symbol_pairs']
        assert o.parameter_names == testing_set['parameter_names']
        assert set(o.parameters.keys()) == set(testing_set['parameters'].keys())
        for k,v in testing_set['parameters'].items():
            o.parameters[k] = v

    @pytest.mark.parametrize("testing_set",testing_sets,ids=testing_names)
    def test__evaluate(self,testing_set):
        r_max = 11.
        N_r = 500
        r = r_max * np.linspace(1,100,N_r)/100

        from pypospack.potential import MorsePotential
        o = MorsePotential(symbols=testing_set['symbols'])
        o.evaluate(r,testing_set['parameters'])

        assert type(o.potential_evaluations) is OrderedDict
        for pair_key,pair_values in o.potential_evaluations.items():
            assert type(pair_values) is np.ndarray
            assert pair_values.shape == r.shape

if __name__ == "__main__":
    from pypospack.potential import MorsePotential

    symbols = ['Ni']
    symbol_pairs = [['Ni','Ni']]
    parameter_names = ['NiNi_D0', 'NiNi_a', 'NiNi_r0']
    parameters = OrderedDict()
    parameters['NiNi_D0'] = 0.001114
    parameters['NiNi_a'] = 3.429506
    parameters['NiNi_r0'] = 2.6813

    r_max = 11.
    N_r = 500
    r = r_max * np.linspace(1,100,N_r)/100

    morse_cls = MorsePotential(symbols=symbols)
    morse_cls.evaluate(r,parameters)
    print("parameters")
    print(80*'=')
    print(morse_cls.parameters)
    print(80*'=')
    print("potential")
    print(80*'=')
    print(morse_cls.potential)

import pytest

from collections import OrderedDict

from pypospack.potential import determine_symbol_pairs
from pypospack.potential import determine_pair_parameter_names
from pypospack.potential import BornMayerPotential

test_cases = OrderedDict()
test_cases['case1'] = OrderedDict()
test_cases['case1']['symbols'] = ['Ni']
test_cases['case1']['potential_type'] = BornMayerPotential.potential_type

parameterized_types = ['symbols','potential_type']
parameterized_test_cases = []
for k,v in test_cases.items():
    
    parameterized_test_cases.append(
            tuple([v[p] for p in parameterized_types])
            )

def test__static_attributes():
    assert BornMayerPotential.potential_type == 'bornmayer'
    assert BornMayerPotential.pair_potential_parameters == ['phi0','gamma','r0']

@pytest.mark.parametrize(
        ','.join(parameterized_types),
        parameterized_test_cases,
        )
def test____init____(symbols,potential_type):
    o = BornMayerPotential(symbols=symbols)
    assert o.potential_type == BornMayerPotential.potential_type
    #assert o.symbol_pairs == determine_symbol_pairs(symbols=symbols)
    assert o.parameter_names == determine_pair_parameter_names(
            symbols=symbols,
            pair_parameter_names=BornMayerPotential.pair_potential_parameters)

if __name__ == "__main__":
    print('test_cases:')
    print(parameterized_test_cases)
    symbols = ['Ni']
    print(BornMayerPotential.potential_type)
    o = BornMayerPotential(symbols=symbols)
    print(o.potential_type)
    print(o.symbol_pairs)
    print(o.parameter_names)

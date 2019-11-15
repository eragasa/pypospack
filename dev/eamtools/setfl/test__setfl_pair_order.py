import pytest
from collections import OrderedDict

from setfl import get_setfl_pair_order

test_cases = OrderedDict()
test_cases['1sym'] = OrderedDict()
test_cases['1sym']['symbols'] = ['Ni']
test_cases['1sym']['expected_pairs'] = [['Ni', 'Ni']]
test_cases['2sym'] = OrderedDict()
test_cases['2sym']['symbols'] = ['Ni','Al']
test_cases['2sym']['expected_pairs'] = [ ['Ni', 'Ni'],
                                         ['Al', 'Ni'],
                                         ['Al', 'Al'] ]


@pytest.mark.parametrize('symbols,expected_pairs',
                         [ (v['symbols'], v['expected_pairs']) for k,v in test_cases.items() ])
def test__setfl_pair_order(symbols,expected_pairs):
    pairs = get_setfl_pair_order(symbols)

    assert isinstance(pairs,list)
    assert all([isinstance(pair,list) for pair in pairs])
    
    for test_pair, expected_pair in zip(pairs,expected_pairs):
        assert test_pair == expected_pair


import pytest

def test__1_symbol():
    from pypospack.potential import determine_symbol_pairs

    symbols = ['Ni']

    pairs = determine_symbol_pairs(symbols)

    assert type(pairs) is list
    assert len(pairs) == 1
    assert pairs[0] == ['Ni','Ni']

def test__2_symbols():
    from pypospack.potential import determine_symbol_pairs

    symbols = ['Mg','O']
    expected_pairs = [['Mg','Mg'],['Mg','O'],['O','O']]
    pairs = determine_symbol_pairs(symbols)

    assert type(pairs) is list
    assert pairs == expected_pairs
    assert pairs[0] == ['Mg','Mg']
    assert pairs[1] == ['Mg','O']
    assert pairs[2] == ['O','O']
    assert len(pairs) == 3

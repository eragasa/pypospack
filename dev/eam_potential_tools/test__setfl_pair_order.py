import pytest

from dev_SeatonSetflReader import get_setfl_pair_order

@pytest.mark.parametrize("symbols,expected_pairs",
        [
            (
                ['Ni'],
                [['Ni','Ni']]
            ),
            (
                ['Ni','Al'],
                [['Ni','Ni'],['Al','Ni'],['Al','Al']]
            )
        ])
def test__get_setfl_pair_order(symbols,expected_pairs):

    pairs = get_setfl_pair_order(symbols)

    assert isinstance(pairs,list)
    assert all([isinstance(pair,list) for pair in pairs])
    
    expected_pairs = [
                ['Ni','Ni'],
                ['Al','Ni'],
                ['Al','Al']
            ]

    for test_pair, expected_pair in zip(pairs,expected_pairs):
        assert test_pair == expected_pair

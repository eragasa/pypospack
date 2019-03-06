import pytest
from pypospack.potential import determine_symbol_pairs

test_sets = [
        (
            'Ni',
            [['Ni','Ni']]
        ),
        (
            ['Ni'],
            [['Ni','Ni']]
        ),
        (
            ['Ni','Al'],
            [
                ['Ni','Ni'],
                ['Ni','Al'],
                ['Al','Al']
            ]
        )
]
test_ids=['Ni_str','Ni_list','NiAl_list']

@pytest.mark.parametrize('symbols,symbol_pairs',test_sets,ids=test_ids)
def test__determine_symbol_pairs(symbols,symbol_pairs):
    assert symbol_pairs == determine_symbol_pairs(symbols=symbols)

if __name__ == "__main__":
    symbols = ['Ni']
    symbol_pairs = determine_symbol_pairs(symbols=symbols)
    print(symbol_pairs)

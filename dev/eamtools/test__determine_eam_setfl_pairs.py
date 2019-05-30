import pytest
from potentialfitter import determine_eam_setfl_pairs

@pytest.mark.parametrize('symbols,expected_pairs',[
    (
        ['Ni'],
        [('Ni','Ni')]
    ),(
        ['Ni','Al'],
        [('Ni','Ni'),('Al','Ni'),('Al','Al')]
    )])
def test__determine_eam_pairs(symbols,expected_pairs):
     pairs = determine_eam_setfl_pairs(symbols)
     assert isinstance(pairs,list)
     assert all([isinstance(k,tuple) for k in pairs])
     assert pairs==expected_pairs

if __name__ == "__main__":
    symbols = ['Ni']
    pairs = determine_eam_setfl_pairs(symbols)
    print(pairs)

    symbols = ['Ni','Al']
    pairs = determine_eam_setfl_pairs(symbols)
    print(pairs)

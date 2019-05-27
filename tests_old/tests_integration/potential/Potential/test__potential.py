import pytest

def test__import__from_pypospack_potential():
    from pypospack.potential import Potential

def test___init____():
    #<--- test vars
    symbols = ['Ni']

    #<--- test setup
    from pypospack.potential import Potential

    #<---- code being tested
    with pytest.raises(NotImplementedError):
        testpot = Potential(symbols=symbols)


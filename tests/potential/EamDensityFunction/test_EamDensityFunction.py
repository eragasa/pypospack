import pytest

def test__import__from_pypospack_potential():
    from pypospack.potential import EamDensityFunction

def test____init__():
    #<--- variables being tested
    symbols = ['Ni']
    #<--- code to setup test
    from pypospack.potential import EamDensityFunction

    #<--- code being tested
    with pytest.raises(NotImplementedError):
        eamdens = EamDensityFunction(symbols=symbols)

import pytest
from collections import OrderedDict

def test__import__from_pypospack_potential():
    from pypospack.potential import BjsEmbeddingFunction

def test_1sym____init__():
    #<--- variables
    symbols = ['Ni']

    #<--- setup
    from pypospack.potential import BjsEmbeddingFunction

    #<--- code being tested
    embed = BjsEmbeddingFunction(symbols=symbols)

    #<--- expected results
    parameter_names = ['Ni_F0','Ni_gamma','Ni_F1']
    
    #<--- testing expected results
    assert type(embed.symbols) is list
    assert len(embed.symbols) == len(symbols)
    assert embed.symbols == symbols

    assert type(embed.parameter_names) is list
    assert len(embed.parameter_names) == len(parameter_names)
    assert embed.parameter_names == parameter_names

    assert type(embed.parameters) is OrderedDict
    assert len(embed.parameters) is len(embed.parameter_names)
    for p in embed.parameter_names:
        assert p in embed.parameters

    assert embed.embedding is None

def test_2sym__init__():
    #<--- variables
    symbols = ['Ni','Al']

    #<--- setup
    from pypospack.potential import BjsEmbeddingFunction

    #<--- code being tested
    embed = BjsEmbeddingFunction(symbols=symbols)

    #<--- expected results
    parameter_names = ['Ni_F0','Ni_gamma','Ni_F1',
                       'Al_F0','Al_gamma','Al_F1']
    
    #<--- testing expected results
    assert type(embed.symbols) is list
    assert len(embed.symbols) == len(symbols)
    assert embed.symbols == symbols

    assert type(embed.parameter_names) is list
    assert len(embed.parameter_names) == len(parameter_names)
    assert embed.parameter_names == parameter_names

    assert type(embed.parameters) is OrderedDict
    assert len(embed.parameters) is len(embed.parameter_names)
    for p in embed.parameter_names:
        assert p in embed.parameters

    assert embed.embedding is None

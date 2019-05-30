import pytest
from collections import OrderedDict
import numpy as np
import pypospack.potential as potential


def test__import__pypospack_potential():
    from pypospack.potential import TersoffPotential

def test__1element____init__():
    symbols = ['Si']

    try:
        tersoff = potential.TersoffPotential(symbols=symbols)
    except:
        pytest.fail()

    potential_type_str = "tersoff"

    assert type(tersoff.potential_type) is str
    assert tersoff.potential_type == potential_type_str


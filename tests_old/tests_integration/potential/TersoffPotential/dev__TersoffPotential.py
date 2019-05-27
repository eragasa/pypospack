import pytest
from collections import OrderedDict
import numpy as np
import pypospack.potential as potential

symbols = ['Si']
pot = potential.TersoffPotential(symbols=symbols)

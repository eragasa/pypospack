# -*= coding: utf-8 -*-
__author__ = "Eugene J. Ragasa"
__copyright__ = "Copyright (C) 2019"
__license__ = "Simplifed BSD License"
__version__ = "1.0"

# 2/18/2019 - EJR
# this object breaks out the parameteric sampling from the original PyposmatMonteCarloSampler.

from pypospack.pyposmat.engines import PyposmatBaseSampler

class PyposmatKdeSampler(PyposmatBaseSampler): pass

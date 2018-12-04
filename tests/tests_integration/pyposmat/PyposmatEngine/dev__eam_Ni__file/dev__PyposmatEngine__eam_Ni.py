import os
from collections import OrderedDict
from pypospack.pyposmat.data import PyposmatDataFile
from pypospack.pyposmat.engines import PyposmatEngine

engine = PyposmatEngine(
        filename_in = 'pyposmat.config.in',
        filename_out = 'pyposmat.config.out')
engine.configure()

_parameters = Ni_eam_parameters
results = engine.evaluate_parameter_set(parameters=_parameters)

print(results)

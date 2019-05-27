import pytest
import os

from collections import OrderedDict

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
parameter_set = OrderedDict()

def test__enforce_parameter_inequality_constraints():
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    sampler = PyposmatBaseSampler()
    sampler.read_configuration_file(filename=config_fn)
    #sampler.enforce_parameter_equality_constraints(free_parameters)
    #sampler.enforce_parameter_inequality_constraints(parameters)

    # tests haven't been implemented yet
    assert False
    

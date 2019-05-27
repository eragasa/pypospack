import pytest 
import os 

from collections import OrderedDict

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
free_parameters = OrderedDict()
parameter_set = OrderedDict()

def test__generate_free_parameters():
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    sampler = PyposmatBaseSampler()
    sampler.read_configuration_file(filename=config_fn)
    #sampler.enforce_parameter_equality_constraints(free_parameters)

    # the tests for this has not been implemented yet
    assert False

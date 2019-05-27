import pytest
import os

from collections import OrderedDict

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')

def test__generate_free_parameters():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

    sampler = PyposmatBaseSampler()
    sampler.read_configuration_file(filename=config_fn)
    free_parameters = sampler.generate_free_parameters()

    assert type(free_parameters) is OrderedDict
    assert all([p in free_parameters for p in sampler.free_parameter_names])
    assert all([p is not None for p in sampler.free_parameter_names])

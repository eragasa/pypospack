import pytest

import os
from pathlib import Path

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')

def test__MgO_configuration_file_exists():
    assert os.path.isdir(pypospack_root_dir)
    assert os.path.isfile(config_fn)

def test__read_configuration_file__set_at__init__statement():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

    config = PyposmatBaseSampler(
            config_fn = config_fn)

def test__read_configuration_file__from_attributes():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

    sampler = PyposmatBaseSampler()
    sampler.config_fn = config_fn
    sampler.read_configuration_file()
    
    do_attribute_tests(sampler=sampler,config_fn=config_fn)

def do_attribute_tests(sampler,config_fn):
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    from pypospack.pyposmat.data import PyposmatConfigurationFile
    
    # test arguments
    assert type(sampler) is PyposmatBaseSampler
    assert type(config_fn) is str


    config = PyposmatConfigurationFile()
    config.read(filename=config_fn)

    assert sampler.structure_directory == config.structures['structure_directory']
    assert sampler.n_iterations == config.sampling_type['n_iterations']
    assert sampler.parameter_names == config.parameter_names
    assert sampler.qoi_names == config.qoi_names
    assert sampler.error_names == config.error_names
    assert sampler.free_parameter_names == config.free_parameter_names
   
    assert set(sampler.parameter_constraints.keys()) == set(config.sampling_constraints.keys())
    assert all([sampler.parameter_constraints[k] == config.sampling_constraints[k] for k in sampler.parameter_constraints])
    assert sampler.parameter_constraints == config.sampling_constraints
    assert sampler.constrained_parameter_names == \
            [p for p in sampler.parameter_names if p not in sampler.free_parameter_names]

if __name__ == "__main__":
    # run this script for output 
    print("pypospack_root_dir:{}".format(pypospack_root_dir))
    print("config_fn:{},{}".format(config_fn,os.path.isfile(config_fn)))

    sampler = PyposmatBaseSampler()
    sampler.config_fn = config_fn
    sampler.read_configuration_file()

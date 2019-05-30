import pytest 
import os 

from collections import OrderedDict

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
config_fn = os.path.join(pypospack_root_dir,'data/MgO_pareto_data/pyposmat.config.in')
parameter_set = OrderedDict()

mpi_rank = 0
mpi_size = 1
log_fn = 'test.log'

def clean_up(log_fn):
    os.remove(log_fn)

def test__run_simulations():
    from pypospack.pyposmat.engines import PyposmatBaseSampler
    
    if True:
        # this is probably what the execution of the PyposmatBaseSampler class looks like
        sampler = PyposmatBaseSampler(
            config_fn = config_fn,
            data_out_fn = data_out_fn,
            mpi_rank = mpi_rank,
            mpi_size = mpi_size)
        sampler.create_base_directories()
        sampler.read_configuration_file()
        sampler.configuration.structures['structure_directory'] = os.path.join('..',_structure_dir)
        sampler.configure_qoi_manager()
        sampler.configure_task_manager()
        sampler.configure_pyposmat_datafile_out()

    assert False
    # the tests for this has not been implemented yet

    assert os.isfile(log_fn)


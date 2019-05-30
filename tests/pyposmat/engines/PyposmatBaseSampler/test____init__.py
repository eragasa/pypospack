import pytest
import os

from pypospack.pyposmat.data import PyposmatLogFile

def test__check_import_from_package_file():
    from pypospack.pyposmat.engines.base_sampler import PyposmatBaseSampler

def test__check_import_from_package__init__():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

def test__check_import_from_package__init__default_arguments():
    from pypospack.pyposmat.engines import PyposmatBaseSampler

    sampler = PyposmatBaseSampler()
    
    assert sampler.mpi_rank == 0
    assert sampler.mpi_size == 1
    assert sampler.config_fn == 'pyposmat.config.in'
    assert type(sampler.data_in_fn) == type(None)
    assert sampler.data_out_fn == 'pyposmat.results.out'
    assert sampler.bad_parameters_fn == 'pyposmat.bad_parameters.out'
    assert sampler.base_directory == os.getcwd()
    assert type(sampler.obj_log) == PyposmatLogFile

if __name__ == "__main__":

    # this area here for developing new tests
    from pypopsack.pyposmat.engines import PyposmatBaseSampler
    sampler = PyposmatBaseSampler()


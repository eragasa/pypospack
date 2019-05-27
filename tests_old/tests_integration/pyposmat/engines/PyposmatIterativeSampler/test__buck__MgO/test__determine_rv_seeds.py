import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')


def test__determine_rv_seeds__no_args__attribute_not_none():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    pyposmat_app.setup_mpi_environment()
    
    pyposmat_app.rv_seed = 1

    assert type(pyposmat_app.rv_seed) is int
    assert type(pyposmat_app.rv_seeds) is type(None)
    assert pyposmat_app.rv_seed == 1

    pyposmat_app.determine_rv_seeds()

    assert type(pyposmat_app.rv_seed) is int
    assert type(pyposmat_app.rv_seeds) is np.ndarray
    assert pyposmat_app.rv_seed == 1
    assert pyposmat_app.rv_seeds.shape[0] == pyposmat_app.mpi_size
    assert pyposmat_app.rv_seeds.shape[1] == pyposmat_app.n_iterations

def test__determine_rv_seeds__seeds_do_not_change_when_run_again():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    pyposmat_app.setup_mpi_environment()
    pyposmat_app.determine_rv_seeds()
    
    # make sure the seed doesn't change when it's run again
    old_seed = pyposmat_app.rv_seed
    old_seeds = pyposmat_app.rv_seeds
    pyposmat_app.determine_rv_seeds()
    
    assert pyposmat_app.rv_seed == old_seed
    assert np.array_equal(pyposmat_app.rv_seeds, old_seeds)

def dev__determine_rv_seeds():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    pyposmat_app.setup_mpi_environment()
    pyposmat_app.determine_rv_seeds()

    print('pyposmat_app.rv_seed:')
    print('\ttype:{}.'.format(str(type(pyposmat_app.rv_seed))))
    print('pyposmat_app.rv_seeds:')
    print(pyposmat_app.rv_seeds.shape)

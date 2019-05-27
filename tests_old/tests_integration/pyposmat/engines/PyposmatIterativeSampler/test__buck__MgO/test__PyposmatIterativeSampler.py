import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')

def test____init__():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename=config_fn)
    pyposmat_app.data_dir = pyposmat_data_dir

    # mpi should not be confiugred
    assert type(pyposmat_app.mpi_comm) is type(None)
    assert type(pyposmat_app.mpi_rank) is type(None)
    assert type(pyposmat_app.mpi_size) is type(None)
    assert type(pyposmat_app.mpi_nprocs) is type(None)

    # testing other attributes aren't set
    assert type(pyposmat_app.n_iterations) is type(None)
    assert type(pyposmat_app.rv_seed) is type(None)
    assert type(pyposmat_app.rv_seeds) is type(None)
    assert type(pyposmat_app.configuration_filename) is str
    assert type(pyposmat_app.configuration) is type(None)
    assert type(pyposmat_app.mc_sampler) is type(None)
    assert type(pyposmat_app.root_directory) is str
    assert type(pyposmat_app.data_directory) is str
    assert type(pyposmat_app.is_restart) is bool
    assert type(pyposmat_app.start_iteration) is int

    # testing values of attributes
    assert pyposmat_app.is_restart == False
    assert pyposmat_app.data_directory == 'data'
    assert pyposmat_app.root_directory == os.getcwd()
    assert pyposmat_app.start_iteration == 0

if False:
    def test__run_all():
        pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
        pyposmat_app.read_configuration_file()
        pyposmat_app.run_all()

if __name__ == "__main__":
    dev__determine_rv_seeds()


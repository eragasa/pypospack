import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')


def test__setup_mpi_environment():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    pyposmat_app.setup_mpi_environment()

    assert type(pyposmat_app.mpi_comm) is MPI.Intracomm
    assert type(pyposmat_app.mpi_rank) is int
    assert type(pyposmat_app.mpi_size) is int

def dev__setup_mpi_environment():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    pyposmat_app.setup_mpi_environment()

    print('pyposmat_app.mpi_rank:{}'.format(str(pyposmat_app.mpi_rank)))
    print('pyposmat_app.mpi_size:{}'.format(str(pyposmat_app.mpi_size)))

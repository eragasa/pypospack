import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')

def test__initialize_rank_directory():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.setup_mpi_environment()
    pyposmat_app.initialize_rank_directory()

    assert os.path.isabs(pyposmat_app.rank_directory)
    assert pyposmat_app.rank_directory == os.path.abspath('rank_{}'.format(pyposmat_app.mpi_rank))



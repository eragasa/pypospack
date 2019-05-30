import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')

def test__read_configuration_file():
    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()

    assert type(pyposmat_app.configuration_filename) is str
    assert type(pyposmat_app.configuration) is PyposmatConfigurationFile
    assert os.path.isabs(pyposmat_app.configuration_filename)
    assert pyposmat_app.configuration_filename == os.path.abspath(config_fn)


import pytest 

import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

pyposmat_data_dir = 'data'
config_fn = os.path.join(pyposmat_data_dir,'pyposmat.config.in')

def test__initialize_data_directory__with_relative_path_arg():
    test_data_dir_path = 'test_dir'

    pyposmat_app = PyposmatIterativeSampler(configuration_filename = config_fn)
    pyposmat_app.read_configuration_file()
    is_created,s = pyposmat_app.initialize_data_directory(data_directory=test_data_dir_path)

    assert s == os.path.join(os.getcwd(),test_data_dir_path)
    assert pyposmat_app.data_directory == os.path.join(os.getcwd(),test_data_dir_path)
    shutil.rmtree(test_data_dir_path)

def test__initialize_data_directory__with_absolute_path_arg():
    test_data_dir_path = os.path.join(os.getcwd(),'test_dir')

    pyposmat_app = PyposmatIterativeSampler(configuration_filename=config_fn)
    pyposmat_app.read_configuration_file()
    is_created,s = pyposmat_app.initialize_data_directory(data_directory=test_data_dir_path)

    assert s == os.path.join(os.getcwd(),test_data_dir_path)
    assert pyposmat_app.data_directory == test_data_dir_path

    shutil.rmtree(test_data_dir_path)

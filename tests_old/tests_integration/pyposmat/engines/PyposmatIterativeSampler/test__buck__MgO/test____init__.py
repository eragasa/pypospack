import pytest
import numpy as np
import os,shutil
import pypospack.utils
from pypospack.pyposmat.data import PyposmatConfigurationFile
from PyposmatIterativeSampler import PyposmatIterativeSampler

pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()
pypospack_data_directory = "./data_test"

config_fn = os.path.join(
        pypospack_root_directory,'MgO_pareto_data')

config_fn = os.path.join
def cleanup():
    if os.isdir(pypospack_data_directory):
        shutil.rmtree(pypospack.data_directory)

def test__init__():
    pyposmat = PyposmatIterativeSampler(configuration_filename=config_fn)

    
if __name__ == "__main__":
    pyposmat = PyposmatIterativeSampler()

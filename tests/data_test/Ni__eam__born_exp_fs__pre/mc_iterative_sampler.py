import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

from pypospack.pyposmat.data import PyposmatDataFile
import pandas as pd
class Dev__PyposmatIterativeSampler(PyposmatIterativeSampler):
    pass 

    
if __name__ == "__main__":
    pyposmat_filename_in = os.path.join(
            'data__Ni__eam__born_exp_fs_00',
            'pyposmat.config.in')
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in)
    pyposmat_app.data_directory = pyposmat_filename_in
    pyposmat_app.read_configuration_file()
    #pyposmat_app.run_restart()
    #exit()
    pyposmat_app.run_all()

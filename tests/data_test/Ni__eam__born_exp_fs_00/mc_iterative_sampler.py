import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

if __name__ == "__main__":
    import Ni__eam__born_exp_fs_0 as config

    pyposmat_filename_in = 'pyposmat.config.in'
 
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in)
    pyposmat_app.read_configuration_file()
    #pyposmat_app.read_configuration_file(filename=pyposmat_configuration_filename)
    pyposmat_app.run_all()

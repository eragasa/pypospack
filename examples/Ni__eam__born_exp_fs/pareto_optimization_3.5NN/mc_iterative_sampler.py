import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

if __name__ == "__main__":
    pyposmat_data_directory = 'data'
    pyposmat_filename_in = os.path.join(
            pyposmat_data_directory,
            'pyposmat.config.in')
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in)
    pyposmat_app.data_directory = pyposmat_data_directory
    pyposmat_app.read_configuration_file()
    pyposmat_app.run_all()

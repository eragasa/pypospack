import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.engines import PyposmatIterativeSampler

# set to true, if you want to restart the simulation, and recover from the last iteration
# which was completed.  
is_restart = False

if __name__ == "__main__":
    pyposmat_data_dir = 'data' 
    pyposmat_filename_in = os.path.join(
            pyposmat_data_dir,'pyposmat.config.in')
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in,
        is_restart=is_restart
    )
    pyposmat_app.data_dir = pyposmat_data_dir
    pyposmat_app.read_configuration_file()
    pyposmat_app.run_all()

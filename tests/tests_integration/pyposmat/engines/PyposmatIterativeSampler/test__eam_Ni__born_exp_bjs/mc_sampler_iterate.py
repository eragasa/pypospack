import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat2.data import PyposmatConfigurationFile
from pypospack.pyposmat2.engines import PyposmatIterativeSampler

if __name__ == "__main__":
    import Ni__eam__morse_exp_fs_0 as config

    pyposmat_filename_in = 'pyposmat.config.in'
    #------------------------------------------------------------------------------
    # WRITE CONFIGURATION FILE
    #------------------------------------------------------------------------------
    configuration = PyposmatConfigurationFile()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    configuration.write(filename=pyposmat_filename_in)

    #------------------------------------------------------------------------------
    # CHECK TO SEE IF CONFIGURATION FILE IS READABLE
    #------------------------------------------------------------------------------
    configuration.read(filename=pyposmat_filename_in)
    
 
    #------------------------------------------------------------------------------
    # RUN PYPOSMAT 
    #------------------------------------------------------------------------------
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pyposmat_filename_in)
    pyposmat_app.read_configuration_file()
    #pyposmat_app.read_configuration_file(filename=pyposmat_configuration_filename)
    pyposmat_app.run_all()

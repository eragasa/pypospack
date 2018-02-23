import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatIterativeSampler
from pypospack.pyposmat import PyposmatDataFile
if __name__ == "__main__":
    import Ni__eam__morse_exp_fs as config

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
    configuration.write(filename='pypospack.config.in')
    configuration.read(filename='pypospack.config.in')
    
    import argparse
    arg_parser = argparse.ArgumentParser()
    arg_parser.add_argument('-i','--input',type=str,
            help="input file to be processed")
    args = vars(arg_parser.parse_args())

    filename = args['input']

    datafile = PyposmatDataFile(filename=filename)
    datafile.read()



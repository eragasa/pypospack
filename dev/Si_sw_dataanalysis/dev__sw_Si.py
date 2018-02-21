import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import PyposmatDataAnalyzer
from pypospack.pyposmat import PyposmatIterativeSampler

if __name__ == "__main__":
    import Si_sw

#WRITE CONFIGURATION FILE

    Si_sw_configuration = PyposmatConfigurationFile()
    Si_sw_configuration.qois = Si_sw.Si_sw_qoi_db.qois
    Si_sw_configuration.potential = Si_sw.potential
    Si_sw_configuration.structures = Si_sw.Si_sw_structures
    Si_sw_configuration.sampling_type = Si_sw.Si_sw_sampling
    Si_sw_configuration.sampling_distribution = Si_sw.Si_sw_param_dist
    Si_sw_configuration.write(filename = 'pypospack.config.in')
    Si_sw_configuration.read(filename = 'pypospack.config.in')

    
    pypospack_filename_in = 'pypospack.config.in'
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = pypospack_filename_in)
    pyposmat_app.read_configuration_file()
    #pyposmat_app.read_configuration_file(filename=pyposmat_configuration_filename)
    pyposmat_app.run_all()

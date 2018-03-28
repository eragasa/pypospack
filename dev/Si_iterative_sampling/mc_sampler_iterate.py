import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.engines import PyposmatIterativeSampler
if __name__ == "__main__":

    _configuration_filename = 'pyposmat.config.in'
    #--------------------------------------------------------------------------
    # DEFINE CONFIGURATION FILED
    #--------------------------------------------------------------------------
    import Si_sw as config
    configuration = PyposmatConfigurationFile()
    configuration.qois = config.qoi_db.qois
    configuration.qoi_constraints = config.qoi_constraints
    configuration.structures = config.structure_db
    configuration.potential = config.potential_formalism
    configuration.sampling_type = config.sampling
    configuration.sampling_distribution = config.parameter_distribution
    configuration.sampling_constraints = config.parameter_constraints
    #<------------- write configuration file
    configuration.write(filename=_configuration_filename)
    #--------------------------------------------------------------------------
    # TEST CONFIGURATION FILE
    #--------------------------------------------------------------------------
    configuration.read(filename=_configuration_filename)
    
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = _configuration_filename)
    pyposmat_app.read_configuration_file()
    pyposmat_app.run_all()

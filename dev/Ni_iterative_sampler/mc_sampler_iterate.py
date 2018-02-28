import os,shutil,sys
import numpy as np
from mpi4py import MPI
from pypospack.pyposmat.data import PyposmatConfigurationFile
from pypospack.pyposmat.data import PyposmatDataAnalyzer
from pypospack.pyposmat.engines import PyposmatIterativeSampler
if __name__ == "__main__":
    import Ni__eam__morse_exp_universal as Ni_eam

    _configuration_filename = 'pyposmat.config.in'
    #------------------------------------------------------------------------------
    # WRITE CONFIGURATION FILE
    #------------------------------------------------------------------------------
    Ni_eam_configuration = PyposmatConfigurationFile()
    Ni_eam_configuration.qois = Ni_eam.Ni_qoi_db.qois
    Ni_eam_configuration.potential = Ni_eam.Ni_eam_potential_formalism
    Ni_eam_configuration.structures = Ni_eam.Ni_structure_db
    Ni_eam_configuration.sampling_type = Ni_eam.Ni_eam_sampling
    Ni_eam_configuration.sampling_distribution =Ni_eam.Ni_eam_parameter_distribution
    Ni_eam_configuration.write(filename=_configuration_filename)
    Ni_eam_configuration.read(filename=_configuration_filename)
    
    pyposmat_app = PyposmatIterativeSampler(
        configuration_filename = _configuration_filename)
    pyposmat_app.read_configuration_file()
    pyposmat_app.run_all()

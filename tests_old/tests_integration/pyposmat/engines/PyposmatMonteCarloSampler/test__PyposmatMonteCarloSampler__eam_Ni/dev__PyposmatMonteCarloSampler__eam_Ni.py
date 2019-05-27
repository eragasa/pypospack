import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import QoiDatabase
from pypospack.qoi import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  
from pypospack.task.lammps import LammpsSimulationError
from pypospack.task.task_manager import PypospackTaskManagerError
import scipy.stats

import Ni__eam__morse_exp_universal as Ni_eam

filename_out='pypospack.results.out'

#------------------------------------------------------------------------------
# WRITE CONFIGURATION FILE
#------------------------------------------------------------------------------
Ni_eam_configuration = PyposmatConfigurationFile()
Ni_eam_configuration.qois = Ni_eam.Ni_qoi_db.qois
Ni_eam_configuration.potential = Ni_eam.Ni_eam_potential_formalism
Ni_eam_configuration.structures = Ni_eam.Ni_structure_db
Ni_eam_configuration.sampling_type = Ni_eam.Ni_eam_sampling
Ni_eam_configuration.sampling_distribution =Ni_eam.Ni_eam_parameter_distribution
Ni_eam_configuration.write(filename='pypospack.config.in')
Ni_eam_configuration.read(filename='pypospack.config.in')
# <---------------- end make configuration file

filename_in='pypospack.config.in'
filename_out='pypospack.results.out'
engine = PyposmatMonteCarloSampler(
        filename_in=filename_in,
        filename_out=filename_out)

# <---------------- printout for debugging purposes
print('base_directory:{}'.format(engine.base_directory))
print('input_filename:{}'.format(engine.pyposmat_filename_in))
print('output_filename:{}'.format(engine.pyposmat_filename_out))
# <---------------- the steps of engine.configure() tested individually
#                   this is the step which configures the object from the
#                   configuration file
# engine.configure()
engine.create_base_directories()
engine.read_configuration_file()
engine.configure_qoi_manager()
engine.configure_task_manager()
engine.configure_pyposmat_datafile_out()
#pyposmat_datafile_out = PyposmatDataFile(filename_out)


engine.print_structure_database()
engine.print_sampling_configuration()
engine.print_initial_parameter_distribution()

for i in range(engine.n_iterations):
    engine.run_simulations(i_iteration=i,n_samples=1000)

exit()


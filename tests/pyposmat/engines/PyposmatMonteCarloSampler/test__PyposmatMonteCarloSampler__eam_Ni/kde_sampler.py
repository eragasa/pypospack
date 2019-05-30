import copy,yaml,os
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

filename_in = os.path.join('data','pypospack.kde.01.out')
engine.run_simulations(i_iteration=1,n_samples=1000,filename=filename_in)
#engine.run_simulations(i_iteration=i,n_samples=1000)

exit()

# print out qoi performance constraints
#print(80*'-')
#print('{:^80}'.format('QOI PERFORMANCE CONSTRAINTS'))
#print(80*'-')
#for qoi_name, value in Si_sw_qoi_constraints.items():
#    print('{:>20} {:>10}'.format(qoi_name,value))

_structure_directory = engine.configuration.structures['structure_directory']
_n_iterations = engine.configuration.sampling_type['n_iterations']
_n_samples = engine.configuration.sampling_type[1]['n_samples']
_rv_generators = OrderedDict()

_filename_in = os.path.join(
        'data',
        'pypospack.kde.01.out')
_datafile_in = PyposmatDataFile(filename=_filename_in)
_datafile_in.read()
_rv_generator = scipy.stats.gaussian_kde(
        _datafile_in.df[engine.free_parameter_names].values.T)



engine.pyposmat_datafile_out.write_header_section(
        filename=filename_out,
        parameter_names=engine.parameter_names,
        qoi_names=engine.qoi_names,
        error_names=engine.error_names)

import time
time_start = time.time()
_n_errors = 0
_n_samples = engine.configuration.sampling_type[1]['n_samples']
for i_sample in range(_n_samples):
    # generate parameter set
    _parameters = OrderedDict([(p,None) for p in engine.parameter_names])
    _free_parameters = _rv_generator.resample(1)
    for i,v in enumerate(engine.free_parameter_names):
        _parameters[v] = _free_parameters[i]

    _constrained_parameters = []
    for p in engine.parameter_names:
        if p not in engine.free_parameter_names:
            _constrainted_parameters.append(p)

    for p in _constrained_parameters:
        _str_eval = str(_param_dist_def[p][1])
        for fp in free_parameter_names:
            if fp in _str_eval:
                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
        _parameters[p] = eval(_str_eval)

    try:
        _results = engine.evaluate_parameter_set(parameters=_parameters)
    except LammpsSimulationError as e:
        _n_errors += 1
    except PypospackTaskManagerError as e:
        _n_errors += 1
    else:
        engine.pyposmat_datafile_out.write_simulation_results(
                filename=filename_out,
                sim_id=i_sample,
                results=_results)
    finally:
        # print out summaries every 10 solutions
        if (i_sample+1)%10 == 0:
            n_samples_completed = i_sample+1
            time_end = time.time()
            time_total = time_end-time_start
            avg_time = n_samples_completed/time_total
            _str_msg = '{} samples completed in {:.4f}s. Avg_time = {:.4f}. n_errors = {}'.format(
                n_samples_completed,
                time_total,
                avg_time,
                _n_errors)
            print(_str_msg)

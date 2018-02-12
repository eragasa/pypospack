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
engine.pyposmat_datafile_out = PyposmatDataFile(filename_out)

_structure_directory = engine.configuration.structures['structure_directory']
_n_iterations = engine.configuration.sampling_type['n_iterations']
_n_samples = engine.configuration.sampling_type[0]['n_samples']

#param_dist_def = Si_sw_param_dist['parameters']
parameter_names = [p for p in engine.configuration.sampling_distribution]
qoi_names = [k for k in engine.configuration.qois]
error_names = ['{}.err'.format(k) for k in qoi_names]
_param_dist_def = engine.configuration.sampling_distribution
free_parameter_names = [k for k,v in _param_dist_def.items() if v[0] != 'equals']

print(80*'-')
print('{:^80}'.format('STRUCTURE DATABASE'))
print(80*'-')
print('structure_directory:{}'.format(_structure_directory))
print('')
print('{:^20} {:^20}'.format('name','filename'))
print('{} {}'.format(20*'-',20*'-'))
for k,v in engine.structures['structures'].items():
    print('{:20} {:20}'.format(k,v))

print(80*'-')
print('{:^80}'.format('SAMPLING CONFIGURATION'))
print(80*'-')
print('{:^10} {:^10} {:^20}'.format(
    'iteration',
    'n_samples',
    'sampling_type'))
print('{} {} {}'.format(10*'-',10*'-',20*'-'))

for i in range(_n_iterations):
    _n_samples = engine.configuration.sampling_type[i]['n_samples']
    _sample_type = engine.configuration.sampling_type[i]['type']
    print('{:^10} {:^10} {:^20}'.format(i,_n_samples,_sample_type))

print(80*'-')
print('{:80}'.format('INITIAL PARAMETER DISTRIBUTION'))
print(80*'-')
for p in _param_dist_def:
    if p in free_parameter_names:
        str_free = 'free'
        print('{:^20} {:^10} {:^10} {:^10} {:^10}'.format(
            p,
            str_free,
            _param_dist_def[p][0],
            _param_dist_def[p][1]['a'],
            _param_dist_def[p][1]['b']))
    else:
        str_free = 'not_free'
        print('{:^20} {:^10}'.format(p,str_free))

# print out qoi performance constraints
#print(80*'-')
#print('{:^80}'.format('QOI PERFORMANCE CONSTRAINTS'))
#print(80*'-')
#for qoi_name, value in Si_sw_qoi_constraints.items():
#    print('{:>20} {:>10}'.format(qoi_name,value))

_rv_generators = OrderedDict()
for p in free_parameter_names:
    if _param_dist_def[p][0] == 'uniform':
        a = _param_dist_def[p][1]['a']
        b = _param_dist_def[p][1]['b']
        _loc = a
        _scale = b-a
        _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
    else:
        pass

i_iteration = 0
print(80*'-')
print('START SIMULATIONS')
print(80*'-')
print('ITERATION={}'.format(i_iteration))
print(80*'-')

engine.pyposmat_datafile_out.write_header_section(
        filename=filename_out,
        parameter_names=parameter_names,
        qoi_names=qoi_names,
        error_names=error_names)

import time
time_start = time.time()
_n_errors = 0
_n_samples = engine.configuration.sampling_type[0]['n_samples']
for i_sample in range(_n_samples):
    # generate parameter set
    _parameters = OrderedDict([(p,None) for p in parameter_names])
    for p in free_parameter_names:
        _parameters[p] = _rv_generators[p].rvs(size=1)[0]

    _constrained_parameters = [
            p for p in _parameters if p not in free_parameter_names]
    for p in _constrained_parameters:
        _str_eval = str(_param_dist_def[p][1])
        for fp in free_parameter_names:
            if fp in _str_eval:
                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
        _parameters[p] = eval(_str_eval)

    #_parameters = Si.Si_dia_swpizzagalli['parameters']
    try:
        _results = engine.evaluate_parameter_set(parameters=_parameters)
    except LammpsSimulationError as e:
        _n_errors += 1
    except PypospackTaskManagerError as e:
        _n_errors += 1
    else:
        print(_results)
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

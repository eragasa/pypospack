import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import QoiDatabase
from pypospack.qoi import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  

import Si_sw

calc_elastic_properties = False
calc_point_defects = True
# <---------------- making a configuration file
Si_sw_qoi_db = QoiDatabase()
Si_sw_qoi_db.add_qoi(
        qoi_name='Si_dia.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Si')]),
        target=5.431)

# <----------------- ELASTIC PROPERTIES
if calc_elastic_properties:
    Si_sw_qoi_db.add_qoi(
            qoi_name='Si_dia.c11',
            qoi_type='c11',
            structures=OrderedDict([('ideal','Si')]),
            target=151.00)
    Si_sw_qoi_db.add_qoi(
            qoi_name='Si_dia.c12',
            qoi_type='c12',
            structures=OrderedDict([('ideal','Si')]),
            target=75.00)
    Si_sw_qoi_db.add_qoi(
            qoi_name='Si_dia.c44',
            qoi_type='c44',
            structures=OrderedDict([('ideal','Si')]),
            target=56.00)
    Si_sw_qoi_db.add_qoi(
            qoi_name='Si_dia.B',
            qoi_type='bulk_modulus',
            structures=OrderedDict([('ideal','Si')]),
            target=100.00)

# <---------------- define potential formalism
Si_sw_potential = OrderedDict()
Si_sw_potential['potential_type'] = 'stillingerweber'
Si_sw_potential['symbols'] = ['Si']
Si_sw_potential['cutoff_global'] = 10.0
# <---------------- Define Sampling Requirements
Si_sw_param_dist = OrderedDict()
Si_sw_param_dist['mc_sampling'] = OrderedDict()
Si_sw_param_dist['mc_sampling']['seed'] = 0
Si_sw_param_dist['mc_sampling']['n_iterations'] = 10

n_iterations = Si_sw_param_dist['mc_sampling']['n_iterations']
n_samples_per_iteration = 100
for i in range(n_iterations):
    Si_sw_param_dist['mc_sampling'][i] = OrderedDict()
    Si_sw_param_dist['mc_sampling'][i]['type'] = 'kde'
    Si_sw_param_dist['mc_sampling'][i]['n_samples'] = n_samples_per_iteration
    Si_sw_param_dist['mc_sampling'][0]['type'] = 'parametric'
#<----------------- determine parameters
Si_sw_param_dist['parameters'] = OrderedDict()
#<----------------- free parameters
# For uniform distributions, 
#     a = is the low of the rnage, 
#     b = is the high of the
Si_sw_param_dist['parameters']['SiSiSi_epsilon']   = ['uniform',{'a': 2.1, 'b':2.2}]
Si_sw_param_dist['parameters']['SiSiSi_sigma'] = ['uniform',{'a': 1.0, 'b':3.0}] 
Si_sw_param_dist['parameters']['SiSiSi_a']    = ['uniform',{'a': 1.5, 'b':2.0}]
Si_sw_param_dist['parameters']['SiSiSi_lambda']    = ['uniform',{'a': 20.0, 'b':32}]
Si_sw_param_dist['parameters']['SiSiSi_gamma']    = ['uniform',{'a': 1.0, 'b':2.0}]
Si_sw_param_dist['parameters']['SiSiSi_costheta0']    = ['equals',-0.333333333333]
Si_sw_param_dist['parameters']['SiSiSi_A']    = ['uniform',{'a': 6.0, 'b':20.0}]
Si_sw_param_dist['parameters']['SiSiSi_B']    = ['uniform',{'a': 0.5, 'b':1.0}]
Si_sw_param_dist['parameters']['SiSiSi_p']    = ['equals',4.0]
Si_sw_param_dist['parameters']['SiSiSi_q']    = ['equals',0.0]
Si_sw_param_dist['parameters']['SiSiSi_tol']    = ['equals',0.0]

#<----------------- constrained parameters
#<----------------- parameter constriants
#MgO_parameter_constraints = OrderedDict()
#MgO_parameter_constraints['chrgMg_gt_0'] = ['chrg_Mg > 0']
#MgO_parameter_constraints['chrgO_lt_0'] = ['chrg_O < 0']
#MgO_parameter_constraints['MgMg_A_gt_0']  = ['MgMg_A > 0']
#MgO_parameter_constraints['MgMg_rho_gt_0']  = ['MgMg_rho > 0']
#MgO_parameter_constraints['MgMg_C_gt_0']  = ['MgMg_C > 0']
#MgO_parameter_constraints['MgO_A_gt_0']  = ['MgO_A > 0']
#MgO_parameter_constraints['MgO_rho_gt_0']  = ['MgO_rho > 0']
#MgO_parameter_constraints['MgO_C_gt_0']  = ['MgO_C > 0']
#MgO_parameter_constraints['OO_A_gt_0']  = ['OO_A > 0']
#MgO_parameter_constraints['OO_rho_gt_0']  = ['OO_rho > 0']
#MgO_parameter_constraints['OO_C_gt_0']  = ['OO_C > 0']
#<----------------- qoi performance constraints
Si_sw_qoi_constraints = OrderedDict()

# define performance constraints as 20 of the qoi target value
for qoi_name, qoi_info in Si_sw_qoi_db.qois.items():
    Si_sw_qoi_constraints[qoi_name] = qoi_info['target'] * 0.20

# print out qoi performance constraints
print(80*'-')
print('{:^80}'.format('QOI PERFORMANCE CONSTRAINTS'))
print(80*'-')

for qoi_name, value in Si_sw_qoi_constraints.items():
    print('{:>20} {:>10}'.format(qoi_name,value))


Si_sw_structures = OrderedDict()
Si_sw_structures['structure_directory'] = 'test__PyposmatMonteCarloSampler'
Si_sw_structures['structures'] = OrderedDict()
Si_sw_structures['structures']['Si'] = 'Si_dia_unit.vasp'
Si_sw_configuration = PyposmatConfigurationFile()
Si_sw_configuration.qois = Si_sw_qoi_db.qois
Si_sw_configuration.potential = Si_sw_potential
Si_sw_configuration.structures = Si_sw_structures
assert isinstance(Si_sw_configuration.configuration,OrderedDict)
Si_sw_configuration.write(filename='pypospack.config.in')
Si_sw_configuration.read(filename='pypospack.config.in')
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

n_iterations = Si_sw_param_dist['mc_sampling']['n_iterations']
n_samples = Si_sw_param_dist['mc_sampling'][0]['n_samples']
param_dist_def = Si_sw_param_dist['parameters']

parameter_names = [p for p in param_dist_def]
free_parameter_names = [k for k,v in param_dist_def.items() if v[0] != 'equals']
for p in param_dist_def:
    if p in free_parameter_names:
        str_free = 'free'
        print('{:^10} {:^10} {:^10} {:^10} {:^10}'.format(
            p,
            str_free,
            param_dist_def[p][0],
            param_dist_def[p][1]['a'],
            param_dist_def[p][1]['b']))
    else:
        str_free = 'not_free'
        print('{:^10} {:^10}'.format(p,str_free))
import scipy.stats


_rv_generators = OrderedDict()
for p in free_parameter_names:
    if param_dist_def[p][0] == 'uniform':
        a = param_dist_def[p][1]['a']
        b = param_dist_def[p][1]['b']
        _loc = a
        _scale = b-a
        _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
    else:
        pass

"""
----Proposed Change Outline----

# _parameters could be a property with getter that handles a 'kde', 'uniform' or other flag
@property
def _parameters():
    return None

@_parameters.getter
def _parameters_getter():
    _rv_generators = OrderedDict()
    for p in free_parameter_names:
        if param_dist_def[p][0] == 'uniform':
            a = param_dist_def[p][1]['a']
            b = param_dist_def[p][1]['b']
            _loc = a
            _scale = b-a
            _rv_generators[p] = scipy.stats.uniform(loc=_loc,scale=_scale)
        elif param_dist_def[p][0] == 'kde':
            # load data_file 
            # 
            obj_kd = scipy.stats.gaussian_kde()
            obj_kde = obj_kd(data_file)
            for i_sample in range(n_samples):
                # generate parameter set
                _parameters = OrderedDict([(p,None) for p in parameter_names])
                for p in free_parameter_names:
                    _rv_generators[p] = obj_kd.resample(size=1)
        else:
            pass

    for i_sample in range(n_samples):
        # generate parameter set
        _parameters = OrderedDict([(p,None) for p in parameter_names])
        for p in free_parameter_names:
            try:
                # the 'uniform' case
                _parameters[p] = _rv_generators[p].rvs(size=1)[0]
            except (some error about .rvs() not being a method of _rv_generators[p]):
                # the 'kde' case
                _parameters[p] = _rv_generators[p]
            else:
                # some error
                _parameters[p] = None

    return _parameters
"""

for i_sample in range(n_samples):
    # generate parameter set
    _parameters = OrderedDict([(p,None) for p in parameter_names])
    for p in free_parameter_names:
        _parameters[p] = _rv_generators[p].rvs(size=1)[0]

    _constrained_parameters = [
            p for p in _parameters if p not in free_parameter_names]
    for p in _constrained_parameters:
        _str_eval = str(param_dist_def[p][1])
        for fp in free_parameter_names:
            if fp in _str_eval:
                _str_eval = _str_eval.replace(fp,str(_parameters[fp]))
        _parameters[p] = eval(_str_eval)

    #_parameters = Si.Si_dia_swpizzagalli['parameters']
    _results = engine.evaluate_parameter_set(parameters=_parameters)
    
    _strout = str(i_sample) + ","\
            + ",".join([str(v) for k,v in _results['parameters'].items()]) + ","\
            + ",".join([str(v) for k,v in _results['qois'].items()]) + ","\
            + ",".join([str(v) for k,v in _results['errors'].items()])
    #print(_strout)
    print(i_sample)
    print(_results)
    print(_results['parameters'])
    print(_results['qois'])
    print(_results['errors'])
    print(_results['parameters']['SiSiSi_a'])

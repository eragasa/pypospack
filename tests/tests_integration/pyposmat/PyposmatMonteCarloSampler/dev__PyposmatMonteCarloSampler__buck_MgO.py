import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
#from pypospack.pyposmat import QoiDatabase
from pypospack.qoi import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  

import MgO

# <---------------- making a configuration file
MgO_qoi_db = QoiDatabase()
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.c11',
#        qoi_type='c11',
#        structures=OrderedDict([('ideal','MgO_NaCl')]),
#        target=277.00)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.c12',
#        qoi_type='c12',
#        structures=OrderedDict([('ideal','MgO_NaCl')]),
#        target=91.67)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.c44',
#        qoi_type='c44',
#        structures=OrderedDict([('ideal','MgO_NaCl')]),
#        target=144.01)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.B',
#        qoi_type='bulk_modulus',
#        structures=OrderedDict([('ideal','MgO_NaCl')]),
#        target=153.45)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.G',
#        qoi_type='shear_modulus',
#        structures=OrderedDict([('ideal','MgO_NaCl')]),
#        target=92.66)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.fr_a',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#            ('defect','MgO_NaCl_fr_a'),
#            ('ideal','MgO_NaCl')]),
#        target=10.978)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.fr_c',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#            ('defect','MgO_NaCl_fr_c'),
#            ('ideal','MgO_NaCl')]),
#        target=8.986)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.sch',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#            ('defect','MgO_NaCl_sch'),
#            ('ideal','MgO_NaCl')]),
#        target=5.067)
#MgO_qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.001s',
#        qoi_type='surface_energy',
#        structures=OrderedDict([
#            ('slab','MgO_NaCl_001s'),
#            ('ideal','MgO_NaCl')]),
#        target=0.05595)

# <---------------- define potential formalism
MgO_potential = OrderedDict()
MgO_potential['potential_type'] = 'buckingham'
MgO_potential['symbols'] = ['Mg','O']
MgO_potential['cutoff_global'] = 10.0
# <---------------- Define Sampling Requirements
MgO_param_dist = OrderedDict()
MgO_param_dist['mc_sampling'] = OrderedDict()
MgO_param_dist['mc_sampling']['seed'] = 0
MgO_param_dist['mc_sampling']['n_iterations'] = 10

n_iterations = MgO_param_dist['mc_sampling']['n_iterations']
n_samples_per_iteration = 100
for i in range(n_iterations):
    MgO_param_dist['mc_sampling'][i] = OrderedDict()
    MgO_param_dist['mc_sampling'][i]['type'] = 'kde'
    MgO_param_dist['mc_sampling'][i]['n_samples'] = n_samples_per_iteration
MgO_param_dist['mc_sampling'][0]['type'] = 'parametric'
#<----------------- determine parameters
MgO_param_dist['parameters'] = OrderedDict()
#<----------------- free parameters
# For uniform distributions, 
#     a = is the low of the rnage, 
#     b = is the high of the
MgO_param_dist['parameters']['chrg_Mg'] = ['uniform',{'a':+1.5,  'b':+2.5}]
MgO_param_dist['parameters']['chrg_O']   = ['equals','-chrg_Mg']
MgO_param_dist['parameters']['MgMg_A']   = ['equals',0.000]
MgO_param_dist['parameters']['MgMg_rho'] = ['equals',0.500] 
MgO_param_dist['parameters']['MgMg_C']    = ['equals',0.000]
MgO_param_dist['parameters']['MgO_A']   = ['uniform',{'a':800.00,'b':1300.00}]
MgO_param_dist['parameters']['MgO_rho'] = ['uniform',{'a':0.2900,'b':0.3300}]
MgO_param_dist['parameters']['MgO_C']    = ['equals',0.000]
MgO_param_dist['parameters']['OO_A']    = ['uniform',{'a':500.00,'b':25000.00}]
MgO_param_dist['parameters']['OO_rho']  = ['uniform',{'a':0.1000,'b':0.4000}]
MgO_param_dist['parameters']['OO_C']    = ['uniform',{'a':25.00, 'b':77.00}]
#<----------------- constrained parameters
#<----------------- parameter constriants
MgO_parameter_constraints = OrderedDict()
MgO_parameter_constraints['chrgMg_gt_0'] = ['chrg_Mg > 0']
MgO_parameter_constraints['chrgO_lt_0'] = ['chrg_O < 0']
MgO_parameter_constraints['MgMg_A_gt_0']  = ['MgMg_A > 0']
MgO_parameter_constraints['MgMg_rho_gt_0']  = ['MgMg_rho > 0']
MgO_parameter_constraints['MgMg_C_gt_0']  = ['MgMg_C > 0']
MgO_parameter_constraints['MgO_A_gt_0']  = ['MgO_A > 0']
MgO_parameter_constraints['MgO_rho_gt_0']  = ['MgO_rho > 0']
MgO_parameter_constraints['MgO_C_gt_0']  = ['MgO_C > 0']
MgO_parameter_constraints['OO_A_gt_0']  = ['OO_A > 0']
MgO_parameter_constraints['OO_rho_gt_0']  = ['OO_rho > 0']
MgO_parameter_constraints['OO_C_gt_0']  = ['OO_C > 0']
#<----------------- qoi performance constraints
MgO_qoi_constraints = OrderedDict()

# define performance constraints as 20 of the qoi target value
for qoi_name, qoi_info in MgO_qoi_db.qois.items():
    MgO_qoi_constraints[qoi_name] = qoi_info['target'] * 0.20

# print out qoi performance constraints
print(80*'-')
print('{:^80}'.format('QOI PERFORMANCE CONSTRAINTS'))
print(80*'-')

for qoi_name, value in MgO_qoi_constraints.items():
    print('{:>20} {:>10}'.format(qoi_name,value))


MgO_structures = OrderedDict()
MgO_structures['structure_directory'] = 'test__PyposmatMonteCarloSampler'
MgO_structures['structures'] = OrderedDict()
MgO_structures['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
MgO_configuration = PyposmatConfigurationFile()
MgO_configuration.qois = MgO_qoi_db.qois
MgO_configuration.potential = MgO_potential
MgO_configuration.structures = MgO_structures
assert isinstance(MgO_configuration.configuration,OrderedDict)
MgO_configuration.write(filename='pypospack.config.in')
MgO_configuration.read(filename='pypospack.config.in')
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

n_iterations = MgO_param_dist['mc_sampling']['n_iterations']
n_samples = MgO_param_dist['mc_sampling'][0]['n_samples']
param_dist_def = MgO_param_dist['parameters']

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

    #_parameters = MgO.MgO_LewisCatlow['parameters']
    _results = engine.evaluate_parameter_set(parameters=_parameters)
    print(i_sample,_results)

    

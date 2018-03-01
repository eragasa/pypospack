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
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
Si_sw_potential = OrderedDict()
Si_sw_potential['potential_type'] = 'stillingerweber'
Si_sw_potential['symbols'] = ['Si']
Si_sw_potential['cutoff_global'] = 10.0

Si_sw_qoi_db = QoiDatabase()
Si_sw_qoi_db.add_qoi(
        qoi_name='Si_dia.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Si')]),
        target=5.431)
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
Si_sw_qoi_db.qoi_qoi(
        qoi_name='Si_dia.vac',
        qoi_type='defect_calculation',
        structures=OrderedDict(
            ['defect','Si_vac'],
            ['ideal','Si']
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 10
sampling['mc_seed'] = 0
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(sampling['n_iterations']):
    sampling[i] = OrderedDict()
    sampling[i]['type'] = 'kde'
    sampling[i]['n_samples'] = 10000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
sampling[0]['type'] = 'parametric'
#<----------------- determine parameters

#<----------------- free parameters
# For uniform distributions, 
#     a = is the low of the rnage, 
#     b = is the high of the
Si_sw_param_dist = OrderedDict()
Si_sw_param_dist['SiSiSi_epsilon'] = ['uniform',{'a': 2.1, 'b':2.2}]
Si_sw_param_dist['SiSiSi_sigma'] = ['uniform',{'a': 1.0, 'b':3.0}] 
Si_sw_param_dist['SiSiSi_a'] = ['uniform',{'a': 1.5, 'b':2.0}]
Si_sw_param_dist['SiSiSi_lambda'] = ['uniform',{'a': 20.0, 'b':32}]
Si_sw_param_dist['SiSiSi_gamma'] = ['uniform',{'a': 1.0, 'b':2.0}]
Si_sw_param_dist['SiSiSi_costheta0'] = ['equals',-1/3.]
Si_sw_param_dist['SiSiSi_A'] = ['uniform',{'a': 6.0, 'b':20.0}]
Si_sw_param_dist['SiSiSi_B'] = ['uniform',{'a': 0.5, 'b':1.0}]
Si_sw_param_dist['SiSiSi_p'] = ['equals',4.0]
Si_sw_param_dist['SiSiSi_q'] = ['equals',0.0]
Si_sw_param_dist['SiSiSi_tol'] = ['equals',0.0]

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
Si_sw_structures['structures']['Si_vac'] = 'Si_vac_333.vasp'


#------------------------------------------------------------------------------
# WRITE CONFIGURATION FILE
#------------------------------------------------------------------------------
Si_sw_configuration = PyposmatConfigurationFile()
Si_sw_configuration.qois = Si_sw_qoi_db.qois
Si_sw_configuration.potential = Si_sw_potential
Si_sw_configuration.structures = Si_sw_structures
Si_sw_configuration.sampling_type  = sampling
Si_sw_configuration.sampling_distribution = Si_sw_param_dist
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
    n_samples = 10000
    engine.run_simulations(i,n_samples)

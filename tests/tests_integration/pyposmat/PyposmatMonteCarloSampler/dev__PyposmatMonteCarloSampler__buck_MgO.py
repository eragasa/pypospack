import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
from pypospack.pyposmat import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  

import MgO

# <---------------- making a configuration file
MgO_qoi_db = QoiDatabase()
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=277.00)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=91.67)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=144.01)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=153.45)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=92.66)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_a',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_a'),
            ('ideal','MgO_NaCl')]),
        target=10.978)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_c',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_c'),
            ('ideal','MgO_NaCl')]),
        target=8.986)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.sch',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_sch'),
            ('ideal','MgO_NaCl')]),
        target=5.067)
MgO_qoi_db.add_qoi(
        qoi_name='MgO_NaCl.001s',
        qoi_type='surface_energy',
        structures=OrderedDict([
            ('slab','MgO_NaCl_001s'),
            ('ideal','MgO_NaCl')]),
        target=0.05595)

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
    MgO_param_dist['mc_sampling'][i] = OrdereDict()
    MgO_param_dist['mc_sampling'][i]['type'] = 'kde'
    MgO_param_dist['mc_sampling'][i]['n_samples'] = n_samples_per_iteration
MgO_param_dist['mc_sampling'][0]['type'] = 'parametric'
#<----------------- determine parameters
MgO_param_dist['parameters'] = OrderedDict()
#<----------------- free parameters
MgO_param_dist['parameters']['chrg_Mg'] = ['uniform',{'a':+1.5,  'b':+2.5}]
MgO_param_dist['parameters']['MgO_A']   = ['uniform',{'a':800.00,'b':1300.00}]
MgO_param_dist['parameters']['MgO_rho'] = ['uniform',{'a':0.2900,'b':0.3300}]
MgO_param_dist['parameters']['OO_A']    = ['uniform',{'a':500.00,'b':25000.00}]
MgO_param_dist['parameters']['OO_rho']  = ['uniform',{'a':0.1000,'b':0.4000}]
MgO_param_dist['parameters']['OO_C']    = ['uniform',{'a':25.00, 'b':77.00}]
#<----------------- constrained parameters
MgO_param_dist['parameters']['chrg_O']   = ['equals','-chrg_Mg']
MgO_param_dist['parameters']['MgMg_A']   = ['equals',0.000]
MgO_param_dist['parmaeters']['MgMg_rho'] = ['equals',0.500] 
MgO_param_dist['parameters']['MgO_C']    = ['equals',0.000]
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

MgO_structures = OrderedDict()
MgO_structures['structure_directory'] = 'test__PyposmatEngine'
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

engine = PyposmatEngine(
        filename_in = 'pypospack.config.in',
        filename_out = 'pypospack.config.out')
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
# <
_parameters = MgO.MgO_LewisCatlow['parameters']
engine.evaluate_parameter_set(parameters=_parameters)


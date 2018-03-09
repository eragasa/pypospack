import copy,yaml
from collections import OrderedDict
from pypospack.pyposmat import PyposmatMonteCarloSampler
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
#from pypospack.pyposmat import QoiDatabase
from pypospack.qoi import QoiDatabase
from pypospack.io.filesystem import OrderedDictYAMLLoader  

import  buck_MgO as MgO

# <---------------- making a configuration file
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=277.00)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=91.67)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=144.01)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=153.45)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=92.66)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_a',
        qoi_type='E_formation_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_a'),
            ('ideal','MgO_NaCl')]),
        target=10.978)
#qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.fr_c',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#            ('defect','MgO_NaCl_fr_c'),
#            ('ideal','MgO_NaCl')]),
#        target=8.986)
#qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.sch',
#        qoi_type='point_defect',
#        structures=OrderedDict([
#            ('defect','MgO_NaCl_sch'),
#            ('ideal','MgO_NaCl')]),
#        target=5.067)
#qoi_db.add_qoi(
#        qoi_name='MgO_NaCl.001s',
#        qoi_type='surface_energy',
#        structures=OrderedDict([
#            ('slab','MgO_NaCl_001s'),
#            ('ideal','MgO_NaCl')]),
#        target=0.05595)

# <---------------- define potential formalism
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'buckingham'
potential_formalism['symbols'] = ['Mg','O']
potential_formalism['cutoff_global'] = 10.0
# <---------------- SAMPLING CONFIGURATION
parameter_sampling = OrderedDict()
parameter_sampling['n_iterations'] = 10
parameter_sampling['mc_seed'] = 0
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(parameter_sampling['n_iterations']):
    parameter_sampling[i] = OrderedDict()
    parameter_sampling[i]['type'] = 'kde'
    parameter_sampling[i]['n_samples'] = 10000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
parameter_sampling[0]['type'] = 'parametric'
#<----------------- determine parameters
parameter_distribution = OrderedDict()
parameter_distribution['chrg_Mg'] = ['uniform',{'a':+1.5,  'b':+2.5}]
parameter_distribution['chrg_O'] = ['equals','-chrg_Mg']
parameter_distribution['MgMg_A'] = ['equals',0.000]
parameter_distribution['MgMg_rho'] = ['equals',0.500] 
parameter_distribution['MgMg_C'] = ['equals',0.000]
parameter_distribution['MgO_A'] = ['uniform',{'a':800.00,'b':1300.00}]
parameter_distribution['MgO_rho'] = ['uniform',{'a':0.2900,'b':0.3300}]
parameter_distribution['MgO_C'] = ['equals',0.000]
parameter_distribution['OO_A'] = ['uniform',{'a':500.00,'b':25000.00}]
parameter_distribution['OO_rho'] = ['uniform',{'a':0.1000,'b':0.4000}]
parameter_distribution['OO_C'] = ['uniform',{'a':25.00, 'b':77.00}]
#<----------------- constrained parameters
#<----------------- parameter constriants
parameter_constraints = OrderedDict()
parameter_constraints['chrgMg_gt_0'] = ['chrg_Mg > 0']
parameter_constraints['chrgO_lt_0'] = ['chrg_O < 0']
parameter_constraints['MgMg_A_gt_0']  = ['MgMg_A > 0']
parameter_constraints['MgMg_rho_gt_0']  = ['MgMg_rho > 0']
parameter_constraints['MgMg_C_gt_0']  = ['MgMg_C > 0']
parameter_constraints['MgO_A_gt_0']  = ['MgO_A > 0']
parameter_constraints['MgO_rho_gt_0']  = ['MgO_rho > 0']
parameter_constraints['MgO_C_gt_0']  = ['MgO_C > 0']
parameter_constraints['OO_A_gt_0']  = ['OO_A > 0']
parameter_constraints['OO_rho_gt_0']  = ['OO_rho > 0']
parameter_constraints['OO_C_gt_0']  = ['OO_C > 0']
#<----------------- qoi performance constraints
qoi_constraints = OrderedDict()
for k,v in qoi_db.qois.items():
    qoi_constraints[k] = abs(v['target'] * 0.20)

def print_qoi_performance_constraints(qoi_constraints):
    _str = 80*'-' + "\n"
    _str += '{:^80}'.format('QOI PERFORMANCE CONSTRAINTS') + "\n"
    _str += 80*'-' + "\n"
    for k,v in qoi_constraints.items():
        _str += '{:>20} {:>10}'.format(qoi_name,value)
    print(_str)

structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
structure_db['structures']['MgO_NaCl_fr_a'] = 'MgO_NaCl_333_fr_a.vasp' 
structure_db['structures']['MgO_NaCl_fr_c'] = 'MgO_NaCl_333_fr_c.vasp'
structure_db['structures']['MgO_NaCl_sch'] = 'MgO_NaCl_333_sch.vasp'
structure_db['structures']['MgO_NaCl_001s'] = 'MgO_NaCl_001s.vasp'

#------------------------------------------------------------------------------
# WRITE CONFIGURATION FILE
#------------------------------------------------------------------------------
configuration = PyposmatConfigurationFile()
configuration.qois = qoi_db.qois
configuration.potential = potential_formalism
configuration.structures = structure_db
configuration.sampling_type = parameter_sampling
configuration.sampling_distribution = parameter_distribution
configuration.parameter_constraints = parameter_constraints
configuration.qoi_constraints = qoi_constraints
configuration.write(filename='pypospack.config.in')
configuration.read(filename='pypospack.config.in')
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


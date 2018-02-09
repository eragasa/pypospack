from collections import OrderedDict
from pypospack.qoi import QoiDatabase

Si_sw_pizzagalli = OrderedDict()
Si_sw_pizzagalli['potential_type'] = 'stillingerweber'
Si_sw_pizzagalli['symbols'] = ['Si']
Si_sw_pizzagalli['parameters']= OrderedDict()

#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------

# <---------------- SAMPLING CONFIGURATION
Si_sw_sampling = OrderedDict()
Si_sw_sampling['n_iterations'] = 10
for i in range(Si_sw_sampling['n_iterations']):
    Si_sw_sampling[i] = OrderedDict()
    Si_sw_sampling[i]['type'] = 'kde'
    Si_sw_sampling[i]['n_samples'] = 10000
Si_sw_sampling[0]['type'] = 'parametric'

# <---------------- POTENTIAL DEFINITION
Si_sw_potential = OrderedDict()
Si_sw_potential['potential_type'] = 'stillingerweber'
Si_sw_potential['symbols'] = ['Si']

# <---------------- INITIAL PARAMETER DEFINITION
Si_sw_parameter_distribution = OrderedDict()
Si_sw_parameter_distribution = OrderedDict()
Si_sw_parameter_distribution['SiSiSi_epsilon'] = ['uniform',{'a': 2.1, 'b':2.2}]
Si_sw_parameter_distribution['SiSiSi_sigma'] = ['uniform',{'a': 1.0, 'b':3.0}] 
Si_sw_parameter_distribution['SiSiSi_a'] = ['uniform',{'a': 1.5, 'b':2.0}]
Si_sw_parameter_distribution['SiSiSi_lambda'] = ['uniform',{'a': 20.0, 'b':32}]
Si_sw_parameter_distribution['SiSiSi_gamma'] = ['uniform',{'a': 1.0, 'b':2.0}]
Si_sw_parameter_distribution['SiSiSi_costheta0'] = ['equals',-0.333333333333]
Si_sw_parameter_distribution['SiSiSi_A'] = ['uniform',{'a': 6.0, 'b':20.0}]
Si_sw_parameter_distribution['SiSiSi_B'] = ['uniform',{'a': 0.5, 'b':1.0}]
Si_sw_parameter_distribution['SiSiSi_p'] = ['equals',4.0]
Si_sw_parameter_distribution['SiSiSi_q'] = ['equals',0.0]
Si_sw_parameter_distribution['SiSiSi_tol'] = ['equals',0.0]

Si_sw_parameter_constraints = OrderedDict()
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

# <---------------- STRUCTURE DATABASE
Si_sw_structures = OrderedDict()
Si_sw_structures['structure_directory'] = 'test__PyposmatMonteCarloSampler'
Si_sw_structures['structures'] = OrderedDict()
Si_sw_structures['structures']['Si_dia'] = 'Si_dia_unit.vasp'

#<00--------------- QOI DATABASE
#Si reference database
# LDA-DFT reference data, originally from:
# Pizzagalli et al.  Phil. Mag A (2008) A 83 1191
# Taken from table
# Pizzagalli et al.  J Phys. Condens. Matter (2013) 055801
Si_sw_qoi_db = QoiDatabase()
# <----------------- STRUCTURAL PROPERTIES
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.Ecoh',
    qoi_type='Ecoh_min_all',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=-4.63)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.a0',
    qoi_type='a11_min_all',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=5.43)
# <----------------- ELASTIC PROPERTIES
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c11',
    qoi_type='c11',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=166.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c12',
    qoi_type='c12',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=64.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c44',
    qoi_type='c44',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=80.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.B',
    qoi_type='bulk_modulus',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=99.0)

# <---------------- QOI PERFORMANCE CONSTRAINTS
# define performance constraints as 20 of the qoi target value
Si_sw_qoi_constraints = OrderedDict()
for qoi_name, qoi_info in Si_sw_qoi_db.qois.items():
    Si_sw_qoi_constraints[qoi_name] = qoi_info['target'] * 0.20

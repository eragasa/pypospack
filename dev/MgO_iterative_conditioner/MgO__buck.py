from collections import OrderedDict
from pypospack.qoi import QoiDatabase
MgO_LewisCatlow = OrderedDict()
MgO_LewisCatlow['potential'] = OrderedDict()
MgO_LewisCatlow['potential']['potential_type'] = 'buckingham'
MgO_LewisCatlow['potential']['symbols'] = ['Mg','O']
MgO_LewisCatlow['parameters'] = OrderedDict()
MgO_LewisCatlow['parameters']['chrg_Mg'] = +2.0
MgO_LewisCatlow['parameters']['chrg_O']  = -2.0
MgO_LewisCatlow['parameters']['MgMg_A']   = 0.0
MgO_LewisCatlow['parameters']['MgMg_rho'] = 0.5
MgO_LewisCatlow['parameters']['MgMg_C']   = 0.0
MgO_LewisCatlow['parameters']['MgO_A']    = 821.6
MgO_LewisCatlow['parameters']['MgO_rho']  = 0.3242
MgO_LewisCatlow['parameters']['MgO_C']    = 0.0
MgO_LewisCatlow['parameters']['OO_A']     = 2274.00
MgO_LewisCatlow['parameters']['OO_rho']   = 0.1490
MgO_LewisCatlow['parameters']['OO_C']     = 27.88

#------------------------------------------------------------------------------
# WRITE STRUCTURES FILE
#------------------------------------------------------------------------------
MgO__buck_structures = OrderedDict()
MgO__buck_structures['structure_directory']= 'structure_db'
MgO__buck_structures['structures']=OrderedDict()
MgO__buck_structures['structures']['MgO_NaCl']= 'MgO_NaCl_unit.gga.relax.vasp'

MgO_qoi_db_gga = QoiDatabase()
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=277.00)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=91.67)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=144.01)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=153.45)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=92.66)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.fr_a',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_a'),
            ('ideal','MgO_NaCl')]),
        target=10.978)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.fr_c',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_c'),
            ('ideal','MgO_NaCl')]),
        target=8.986)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.sch',
        qoi_type='point_defect',
        structures=OrderedDict([
            ('defect','MgO_NaCl_sch'),
            ('ideal','MgO_NaCl')]),
        target=5.067)
MgO_qoi_db_gga.add_qoi(
        qoi_name='MgO_NaCl.001s',
        qoi_type='surface_energy',
        structures=OrderedDict([
            ('slab','MgO_NaCl_001s'),
            ('ideal','MgO_NaCl')]),
        target=0.05595)


#------------------------------------------------------------------------------
# WRITE SAMPLING FILE
#------------------------------------------------------------------------------

# <---------------- SAMPLING CONFIGURATION
MgO_buck_sampling = OrderedDict()
MgO_buck_sampling['n_iterations'] = 10
MgO_buck_sampling['mc_seed'] = None
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(MgO_buck_sampling['n_iterations']):
    MgO_buck_sampling[i] = OrderedDict()
    MgO_buck_sampling[i]['type'] = 'kde'
    MgO_buck_sampling[i]['n_samples'] = 1000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
MgO_buck_sampling[0]['type'] = 'parametric'

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

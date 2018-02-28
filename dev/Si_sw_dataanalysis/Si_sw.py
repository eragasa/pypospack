from collections import OrderedDict
from pypospack.qoi import QoiDatabase

#-----------------------------------------------------------------------------
# DEFINE STRUCTURES
#-----------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'Si_structures'
structure_db['structures'] = OrderedDict()
structure_db['structures']['Si_dia'] = 'Si_dia_unit.vasp'

#  in this case 'Si_dia' is the short name which the rest of the file
#  refers to the structure file which is contained in
#  Si_structures/Si_dia_unit.vasp
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'stillingerweber'
potential_formalism['symbols'] = ['Si']
#potential_formalism['cutoff_global'] = 10.0

#-----------------------------------------------------------------------------
# QOI DEFINITIONS
#-----------------------------------------------------------------------------
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='Si_dia.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=5.431)
qoi_db.add_qoi(
        qoi_name='Si_dia.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=75.00)
qoi_db.add_qoi(
        qoi_name='Si_dia.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=56.00)
qoi_db.add_qoi(
        qoi_name='Si_dia.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=100.00)
#<----------------- qoi performance constraints
qoi_constraints = OrderedDict()
for qoi_name, qoi_info in qoi_db.qois.items():
    qoi_constraints[qoi_name] = qoi_info['target'] * 0.20
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 10
sampling['mc_seed'] = None
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
parameter_distribution = OrderedDict()
parameter_distribution['SiSiSi_epsilon'] = ['uniform',{'a': 2.1, 'b':2.2}]
parameter_distribution['SiSiSi_sigma'] = ['uniform',{'a': 1.0, 'b':3.0}] 
parameter_distribution['SiSiSi_a'] = ['uniform',{'a': 1.5, 'b':2.0}]
parameter_distribution['SiSiSi_lambda'] = ['uniform',{'a': 20.0, 'b':32}]
parameter_distribution['SiSiSi_gamma'] = ['uniform',{'a': 1.0, 'b':2.0}]
parameter_distribution['SiSiSi_costheta0'] = ['equals',-1/3.]
parameter_distribution['SiSiSi_A'] = ['uniform',{'a': 6.0, 'b':20.0}]
parameter_distribution['SiSiSi_B'] = ['uniform',{'a': 0.5, 'b':1.0}]
parameter_distribution['SiSiSi_p'] = ['equals',4.0]
parameter_distribution['SiSiSi_q'] = ['equals',0.0]
parameter_distribution['SiSiSi_tol'] = ['equals',0.0]
#<----------------- parameter constraints
parameter_constraints = OrderedDict()



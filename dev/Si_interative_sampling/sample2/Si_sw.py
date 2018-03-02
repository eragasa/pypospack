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
for k,v in qoi_db.qois.items():
    qoi_constraints['{}.err'.format(k)] = 0.20 * abs(qoi_db.qois[k]['target'])
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
parameter_distribution['SiSiSi_epsilon'] = ['uniform',{'a': 2.05, 'b':2.25}]
parameter_distribution['SiSiSi_sigma'] = ['uniform',{'a': 1.8, 'b':2.4}] 
parameter_distribution['SiSiSi_a'] = ['uniform',{'a': 1.6, 'b':2.4}]
parameter_distribution['SiSiSi_lambda'] = ['uniform',{'a': 22.0, 'b':36.0}]
parameter_distribution['SiSiSi_gamma'] = ['uniform',{'a': 0.54, 'b':1.1}]
parameter_distribution['SiSiSi_costheta0'] = ['equals',-1/3.]
parameter_distribution['SiSiSi_A'] = ['uniform',{'a': 5.0, 'b':15.0}]
parameter_distribution['SiSiSi_B'] = ['uniform',{'a': 0.4, 'b':1.2}]
parameter_distribution['SiSiSi_p'] = ['equals',4.0]
parameter_distribution['SiSiSi_q'] = ['equals',0.0]
parameter_distribution['SiSiSi_tol'] = ['equals',0.0]
#<----------------- parameter constraints
parameter_constraints = OrderedDict()



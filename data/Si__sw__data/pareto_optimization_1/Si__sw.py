from collections import OrderedDict
from pypospack.qoi import QoiDatabase

#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 20
sampling['mc_seed'] = None
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(sampling['n_iterations']):
    sampling[i] = OrderedDict()
    sampling[i]['type'] = 'kde'
    sampling[i]['n_samples'] = 10000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
sampling[0]['type'] = 'parametric'

#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'stillingerweber'
potential_formalism['symbols'] = ['Si']
#potential_formalism['cutoff_global'] = 10.0

#-----------------------------------------------------------------------------
# DEFINE INITIAL PARAMETER DISTRIBUTION
#-----------------------------------------------------------------------------
# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
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
parameter_distribution['SiSiSi_p'] = ['uniform',{'a':3.0,'b':4.4}]
parameter_distribution['SiSiSi_q'] = ['uniform',{'a':0.0,'b':1.0}]
parameter_distribution['SiSiSi_tol'] = ['equals',0.0]

#------------------------------------------------------------------------------
# PARAMETER CONSTRAINTS
#------------------------------------------------------------------------------
#<----------------- parameter constraints
parameter_constraints = OrderedDict()
# No parameter constraints

#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['Si_dia'] = 'Si_dia_unit.vasp'
structure_db['structures']['Si_vac'] = 'Si_dia_333_vac.vasp'
#  in this case 'Si_dia' is the short name which the rest of the file
#  refers to the structure file which is contained in
#  Si_structures/Si_dia_unit.vasp

#-----------------------------------------------------------------------------
# QOI DEFINITIONS
#-----------------------------------------------------------------------------
# L. Pizzagalli et al, J. J. Phys.: Condens. Matter 25 (2013) 055801
# doi:10.1088/0953-8984/25/5/055801
# reference values from Table 2 and Table 3
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='Si_dia.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=-4.63)
qoi_db.add_qoi(
        qoi_name='Si_dia.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=5.43)
qoi_db.add_qoi(
        qoi_name="Si_dia.c11",
        qoi_type='c11',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=166.)
qoi_db.add_qoi(
        qoi_name='Si_dia.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=64.00)
qoi_db.add_qoi(
        qoi_name='Si_dia.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=80.00)
qoi_db.add_qoi(
        qoi_name='Si_dia.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','Si_dia')]),
        target=99.00)
qoi_db.add_qoi(
        qoi_name='Si_dia.vac',
        qoi_type='E_formation',
        structures=OrderedDict([
            ('defect','Si_vac'),
            ('ideal','Si_dia')]),
        target=3.6)
#------------------------------------------------------------------------------
# QOI CONSTRAINTS
# QOI constraints are performed in the order they are iterated through in
# in the dictionary.
#
# If you want to implement new constraints, they should be implemented in
# pypospack.pyposmat.data.DataAnalyzer
#     filter_by_qoi_err:
#          key - the qoi_name as in the qoi_db
#          value - the maximum allowable absolute error
#     filter_by_pareto:
#          filters out dominated points if set to True
#------------------------------------------------------------------------------
#<----------------- qoi performance constraints
qoi_constraints = OrderedDict()
qoi_constraints['qoi_constraints'] = OrderedDict()
qoi_constraints['qoi_constraints']['Si_dia.c11'] = ['>',0.0]
qoi_constraints['qoi_constraints']['Si_dia.c12'] = ['>',0.0]
qoi_constraints['qoi_constraints']['Si_dia.c44'] = ['>',0.0]
#for qoi_name, qoi_info in qoi_db.qois.items():
#    qoi_constraints[qoi_name] = abs(qoi_info['target']) * 0.20
qoi_constraints['filter_by_pareto_membership'] = True
qoi_constraints['filter_by_cost_function'] = OrderedDict([
    ('weighting_scheme_type','scale_by_qoi_target'),
    ('loss_function_type','abs_error'),
    ('cost_function_type','weighted_sum'),
    ('pct_to_keep',0.95),
    ('n_potentials_min',500),
    ('n_potentials_max',10000)
])

latex_labels = OrderedDict()
#------------------------------------------------------------------------------
# WRITE CONFIGURATION FILE
# this is currently creating a race condition, where the file is being written
# by multiple ranks to the same location.
#------------------------------------------------------------------------------
if __name__ == '__main__':
    from pypospack.pyposmat.data import PyposmatConfigurationFile
    pyposmat_filename_in = 'pyposmat.config.in'
    configuration = PyposmatConfigurationFile()
    configuration.qois = qoi_db.qois
    configuration.qoi_constraints = qoi_constraints
    configuration.structures = structure_db
    configuration.potential = potential_formalism
    configuration.sampling_type = sampling
    configuration.sampling_distribution = parameter_distribution
    configuration.sampling_constraints = parameter_constraints
    configuration.write(filename=pyposmat_filename_in)

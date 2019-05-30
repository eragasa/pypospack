import os
from collections import OrderedDict
from pypospack.qoi import QoiDatabase

import pypospack.utils

pypospack_root_dir = pypospack.utils.get_pypospack_root_directory()
#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 1
sampling['mc_seed'] = None
sampling[0] = OrderedDict()
sampling[0]['type'] = 'from_file'
sampling[0]['file'] = os.path.join('data','pyposmat.reference.in')
sampling[0]['n_samples'] = 4
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
structure_db['structure_directory'] = os.path.join(
        pypospack_root_dir,'data','Si__structure_db'
)
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

#------------------------------------------------------------------------------
# QOI VALIDATION DATABASE
#------------------------------------------------------------------------------
qoi_validation_db = QoiDatabase()
qoi_validation_db.add_qoi(
        qoi_name='Si_dia.optinal_ph',
        qoi_type='gamma_phonon_4',
        structures=OrderedDict([('ideal','Si_dia_pum')]),
        target=0
        )
qoi_validation_db.add_qoi(
        qoi_name='Si_dia.th_exp',
        qoi_type='thermal_expansion_coefficient',
        structures=OrderedDict([('ideal','MgO_NaCl')])
        target=000
        qoi_options={
            'temperature_min':0,
            'temperature_max':1000,
            'temperature_step':200,
            'time_total':10,
            'time_step':0.001
            'supercell':[10,10,10]
            }
        )
qoi_validation_db.add_qoi(
        qoi_name='Si_dia.rdf',
        qoi_type='rdf',
        structures=OrderedDict([('ideal','Si_dia_unit')])
        target=
        qoi_options={
            'temperature_min':0,
            'temperature_max':3000,
            'temperature_step':200,
            'time_total':10,
            'time_step':0.001,
            'supercell':[10,10,10]
            }
        )
latex_labels = OrderedDict()
#------------------------------------------------------------------------------
# REFERENCE POTENTIALS
#------------------------------------------------------------------------------
reference_potentials=OrderedDict()
reference_potentials['SW'] = OrderedDict()
reference_potentials['SW']['formalism'] = potential_formalism
reference_potentials['SW']['parameters'] = OrderedDict([
    ('SiSiSi_epsilon',2.1686),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',21.0),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])
reference_potentials['VBWM'] = OrderedDict()
reference_potentials['VBWM']['formalism'] = potential_formalism
reference_potentials['VBWM']['parameters'] = OrderedDict([
    ('SiSiSi_epsilon',1.64833),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.5),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])
reference_potentials['PG'] = OrderedDict()
reference_potentials['PG']['formalism'] = potential_formalism
reference_potentials['PG']['parameters'] = OrderedDict([
    ('SiSiSi_epsilon',1.04190),
    ('SiSiSi_sigma',2.128117),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.0),
    ('SiSiSi_gamma',1.10),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',19.0),
    ('SiSiSi_B',0.65),
    ('SiSiSi_p',3.5),
    ('SiSiSi_q',0.5),
    ('SiSiSi_tol',0.0)
])

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

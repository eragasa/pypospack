from collections import OrderedDict
from pypospack.qoi import QoiDatabase

""" Short commment

LONG COMMMENT

"""

__author__ = "Manuel Esparragoza, 2018"
__version__ = 1.0

#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 1
sampling['mc_seed'] = None
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(sampling['n_iterations']):
    sampling[i] = OrderedDict()
    sampling[i]['type'] = 'kde'
    sampling[i]['n_samples'] = 50000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
# ---- Example for parametric sampling
# sampling[0]['type'] = 'parametric'
# ---- Example for file sampling
sampling[0]['type'] = 'from_file'
sampling[0]['file'] = 'data/pyposmat.kde.00.out' 
# ---- Example for KDE sampling
# sampling[0]['type'] = 'kde'
# sampling[0]['file'] = 'data/pyposmat.kde.00.out'
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'buckingham'
potential_formalism['symbols'] = ['Mg','O']

# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
parameter_distribution = OrderedDict()
parameter_distribution['chrg_Mg'] = ['uniform',{'a':+1.5,  'b':+2.5}]
parameter_distribution['chrg_O']   = ['equals','-chrg_Mg']
parameter_distribution['MgMg_A']   = ['equals',0.000]
parameter_distribution['MgMg_rho'] = ['equals',0.500] 
parameter_distribution['MgMg_C']    = ['equals',0.000]
parameter_distribution['MgO_A']   = ['uniform',{'a':800.00,'b':1300.00}]
parameter_distribution['MgO_rho'] = ['uniform',{'a':0.2900,'b':0.3300}]
parameter_distribution['MgO_C']    = ['equals',0.000]
parameter_distribution['OO_A']    = ['uniform',{'a':500.00,'b':25000.00}]
parameter_distribution['OO_rho']  = ['uniform',{'a':0.1000,'b':0.4000}]
parameter_distribution['OO_C']    = ['uniform',{'a':25.00, 'b':77.00}]
#------------------------------------------------------------------------------
# PARAMETER CONSTRAINTS
#------------------------------------------------------------------------------
parameter_constraints = OrderedDict()
parameter_constraints['chrg_Mg > 0'] = 'chrg_Mg > 0.'
parameter_constraints['MgO_A > 0'] = 'MgO_A > 0'
parameter_constraints['MgO_rho > 0'] = 'MgO_rho > 0.'
parameter_constraints['OO_A > 0'] = 'OO_A > 0.'
parameter_constraints['OO_rho > 0'] = 'OO_rho > 0.'
parameter_constraints['OO_C > 0'] = 'OO_C > 0.'
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
structure_db['structures']['MgO_NaCl_fr_a'] = 'MgO_NaCl_333_fr_a.vasp' 
structure_db['structures']['MgO_NaCl_fr_c'] = 'MgO_NaCl_333_fr_c.vasp'
structure_db['structures']['MgO_NaCl_sch'] = 'MgO_NaCl_333_sch.vasp'
structure_db['structures']['MgO_NaCl_001s'] = 'MgO_NaCl_001s.vasp'
#------------------------------------------------------------------------------
# FITTING DATABASE
#------------------------------------------------------------------------------
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.a0',
        qoi_type='a1_min_all',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246)
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
qoi_constraints = OrderedDict()
#qoi_constraints['qoi_constraints']=OrderedDict()
#qoi_constraints['select_pareto_only'] = True
#qoi_constraints['filter_by_percentile'] = [80,'pct']

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

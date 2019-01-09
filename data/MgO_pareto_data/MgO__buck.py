from collections import OrderedDict
from pypospack.qoi import QoiDatabase

# This script does not recreate the simulation results contained in this directory.
# These simulations were done using a previous iteration of this software called pyflamestk
# This script is written so that reporting tools can be used on the previous simulation work
# -- EJR, Aug 2018

import pypospack.utils
_pypospack_root_directory = pypospack.utils.get_pypospack_root_directory()

import os
_pyposmat_structure_directory = os.path.join(_pypospack_root_directory,'data','MgO_structure_db')

#------------------------------------------------------------------------------
# CONFIGURAITON SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 10
sampling['mc_seed'] = None
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(sampling['n_iterations']):
    sampling[i] = OrderedDict()
    sampling[i]['type'] = 'kde'
    sampling[i]['n_samples'] = 1000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
sampling[0]['type'] = 'kde'
sampling[0]['file'] = 'data/pyposmat.kde.0.out'
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'buckingham'
potential_formalism['symbols'] = ['Mg','O']

# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
parameter_distribution = OrderedDict()

# this is setup for the assumption of charge neutrality
parameter_distribution['chrg_Mg'] = ['uniform',{'a':+1.5,'b':+2.5}]
parameter_distribution['chrg_O'] = ['equals','-chrg_Mg']

# Ignoring the MgMg term.
parameter_distribution['MgMg_A'] = ['equals', 0.0]
parameter_distribution['MgMg_rho'] = ['equals', 0.5]
parameter_distribution['MgMg_C'] = ['equals',0.0]

# MgO_C == 0, on the assumption of no van der Waals forces
parameter_distribution['MgO_A'] = ['uniform',{'a':800.0,'b':1300.0}]
parameter_distribution['MgO_rho'] = ['uniform',{'a':0.29,'b':0.33}]
parameter_distribution['MgO_C'] = ['equals', 0.0]

parameter_distribution['OO_A'] = ['uniform',{'a':500,'b':25000.0}]
parameter_distribution['OO_rho'] = ['uniform',{'a':0.1,'b':0.4}]
parameter_distribution['OO_C'] = ['uniform',{'a':25.0,'b':77.0}]

parameter_constraints = OrderedDict()
parameter_constraints['chrg_Mg > 0'] = 'chrg_Mg > 0.'
parameter_constraints['chrg_O < 0'] = 'chrg_O < 0.'
parameter_constraints['MgO_A > 0'] = 'MgMg_A > 0.'
parameter_constraints['MgO_rho > 0'] = 'MgMg_rho > 0.'
parameter_constraints['OO_A > 0'] = 'OO_A > 0.'
parameter_constraints['OO_rho > 0'] = 'OO_rho > 0.'
parameter_constraints['OO_C > 0'] = 'OO_C > 0.'
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = _pyposmat_structure_directory
structure_db['structures'] = OrderedDict()
structure_db['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.vasp'
structure_db['structures']['MgO_NaCl_001s'] = 'MgO_NaCl_001s.vasp'
structure_db['structures']['MgO_NaCl_fr_a'] = 'MgO_NaCl_333_fr_a.vasp'
structure_db['structures']['MgO_NaCl_fr_c'] = 'MgO_NaCl_333_fr_c.vasp'
structure_db['structures']['MgO_NaCl_sch'] = 'MgO_NaCl_333_sch.vasp'
#------------------------------------------------------------------------------
# FITTING DATABASE
#------------------------------------------------------------------------------
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.a0",
        qoi_type="a11_min_all",
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=4.246
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.c11",
        qoi_type="c11",
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=277.0
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.c12",
        qoi_type="c12",
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=91.67
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.c44",
        qoi_type="c44",
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=144.01
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.fr_a",
        qoi_type="E_formation",
        structures=OrderedDict(
            [
                ('defect','MgO_NaCl_fr_a'),
                ('ideal','MgO_NaCl')
            ]
        ),
        target=10.978
)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_c',
        qoi_type="E_formation",
        structures=OrderedDict(
            [
                ('defect','MgO_NaCl_fr_c'),
                ('ideal','MgO_NaCl')
            ]
        ),
        target=8.986
)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.sch',
        qoi_type='E_formation',
        structures=OrderedDict(
            [
                ('defect','MgO_NaCl_sch'),
                ('ideal','MgO_NaCl')
            ]
        ),
        target=5.067
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.B",
        qoi_type="bulk_modulus",
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=153.45
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.G",
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','MgO_NaCl')]),
        target=92.66
)
qoi_db.add_qoi(
        qoi_name="MgO_NaCl.001s",
        qoi_type="E_surface",
        structures=OrderedDict(
            [
                ('slab','MgO_NaCl_001s'),
                ('ideal','MgO_NaCl'),
            ],
        ),
        target=0.05595
)

qoi_constraints = OrderedDict()
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

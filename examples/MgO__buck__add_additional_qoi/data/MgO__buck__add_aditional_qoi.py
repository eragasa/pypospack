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
# sampling[0]['type'] = 'from_file'
# sampling[0]['file'] = 'data/pyposmat.kde.00.out' 
# sample from a kde
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
structure_db['structures']['MgO_NaCl_prim'] = 'MgO_NaCl_prim.gga.relax.vasp'
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
        qoi_type='E_formation',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_a'),
            ('ideal','MgO_NaCl')]),
        target=10.978)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.fr_c',
        qoi_type='E_formation',
        structures=OrderedDict([
            ('defect','MgO_NaCl_fr_c'),
            ('ideal','MgO_NaCl')]),
        target=8.986)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.sch',
        qoi_type='E_formation',
        structures=OrderedDict([
            ('defect','MgO_NaCl_sch'),
            ('ideal','MgO_NaCl')]),
        target=5.067)
qoi_db.add_qoi(
        qoi_name='MgO_NaCl.001s',
        qoi_type='E_surface',
        structures=OrderedDict([
            ('slab','MgO_NaCl_001s'),
            ('ideal','MgO_NaCl')]),
        target=0.05595)
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
qoi_constraints['select_pareto_only'] = True
#qoi_constraints['filter_by_percentile'] = [80,'pct']

#------------------------------------------------------------------------------
# QOI VALIDATION DATABASE
#------------------------------------------------------------------------------
qoi_validation_db = QoiDatabase()
qoi_validation_db.add_qoi(
        qoi_name='MgO_NaCl.optical_ph',
        qoi_type='gamma_phonon_4',
        structures=OrderedDict([
            ('ideal','MgO_NaCl_prim')  # calculate from primitive cell
            ]),
        target=385.
        )
qoi_validation_db.add_qoi(
        qoi_name='MgO_NaCl.th_exp',
        qoi_type='thermal_expansion_coefficient',
        structures=OrderedDict([
            ('ideal','MgO_NaCl')       # calculated from the unit cell
            ]),
        target=10.8e-6,                # per deg C, 
                                       # https://www.crystran.co.uk/optical-materials/magnesium-oxide-mgo
        qoi_options={
            'temperature_min':0,       # in degrees Kelvin
            'temperature_max': 1000,    # in degrees Kelvin
            'temperature_step':200,    # in degrees Kelvin
            'time_total':10,           # 10 picoseconds
            'time_step':0.001,         # 1 femtosecond
            'supercell':[5,5,5]        # 5x5x5 repeat of the unit cell
            }
        )
#------------------------------------------------------------------------------
# REFERENCE POTENTIALS
#------------------------------------------------------------------------------
reference_potentials = OrderedDict()
reference_potentials['LC'] = OrderedDict()
reference_potentials['LC']['formalism'] = potential_formalism
reference_potentials['LC']['parameters'] = OrderedDict([
    ('chrg_Mg',+2.00),('chrg_O',-2.00),
    ('MgMg_A',0.0),('MgMg_rho',0.5),('MgMg_C',0.0),
    ('MgO_A',821.60),('MgO_rho',0.3200),('MgO_C',0.00),
    ('OO_A',22764.00),('OO_rho',0.14900),('OO_C',27.88)
    ])
reference_potentials['LC']['qoi'] = OrderedDict([
    ('MgO_NaCl.a0',4.2107902169152),
    ('MgO_NaCl.c11',307.571810176713),
    ('MgO_NaCl.c12',171.13560278936),
    ('MgO_NaCl.c44',168.168424864137),
    ('MgO_NaCl.B',216.61433858514434),
    ('MgO_NaCl.G',68.21810369367648),
    ('MgO_NaCl.fr_a',9.67958273889598),
    ('MgO_NaCl.fr_c',9.81003415099621),
    ('MgO_NaCl.sch',5.796822583346511),
    ('MgO_NaCl.001s',0.06783775861630922)
    ])
reference_potentials['BG1'] = OrderedDict()
reference_potentials['BG1']['formalism'] = potential_formalism
reference_potentials['BG1']['parameters'] = OrderedDict([
    ('chrg_Mg',+2.00),('chrg_O',-2.00),
    ('MgMg_A',0.0),('MgMg_rho',0.5),('MgMg_C',0.0),
    ('MgO_A',1279.69),('MgO_rho',0.29969),('MgO_C',0.00),
    ('OO_A',9547.96),('OO_rho',0.21916),('OO_C',32.00)
    ])
reference_potentials['BG1']['qoi'] = OrderedDict([
    ('MgO_NaCl.a0',4.20883001371435),
    ('MgO_NaCl.c11',360.10697476092),
    ('MgO_NaCl.c12',162.31431501690),
    ('MgO_NaCl.c44',160.683383081696),
    ('MgO_NaCl.B',228.2452015982429),
    ('MgO_NaCl.G',98.89632987200999),
    ('MgO_NaCl.fr_a',12.428466047278562),
    ('MgO_NaCl.fr_c',11.87773342645869),
    ('MgO_NaCl.sch',7.20095386800221),
    ('MgO_NaCl.001s',0.08064679409333339)
    ])

reference_potentials['BG2'] = OrderedDict()
reference_potentials['BG2']['formalism'] = potential_formalism
reference_potentials['BG2']['parameters'] = OrderedDict([
    ('chrg_Mg',+1.7),('chrg_O',-1.7),
    ('MgMg_A',0.0),('MgMg_rho',0.5),('MgMg_C',0.0),
    ('MgO_A',929.69),('MgO_rho',0.29909),('MgO_C',0.00),
    ('OO_A',4870.0),('OO_rho',0.2670),('OO_C',7.00)
    ])
reference_potentials['BG2']['qoi'] = OrderedDict([
    ('MgO_NaCl.a0',4.222448),
    ('MgO_NaCl.c11',301.31582205825),
    ('MgO_NaCl.c12',150.827961278274),
    ('MgO_NaCl.c44',142.471471981922),
    ('MgO_NaCl.B',200.99058153826635),
    ('MgO_NaCl.G',75.2439303899885),
    ('MgO_NaCl.fr_a',10.43573261594292),
    ('MgO_NaCl.fr_c',8.526618652243087),
    ('MgO_NaCl.sch',5.509124492308274),
    ('MgO_NaCl.001s',0.0692527122209811)
    ])
#------------------------------------------------------------------------------
# LATEX LABELS
#------------------------------------------------------------------------------
latex_labels = OrderedDict()
latex_labels['chrg_Mg'] = OrderedDict([
    ('name',r'Z_{Mg}')])
latex_labels['chrg_O'] = OrderedDict([
    ('name',r'Z_{O}')])

latex_labels['MgMg_A'] = OrderedDict([
    ('name',r'A_{Mg,Mg}')])
latex_labels['MgMg_rho'] = OrderedDict([
    ('name',r'\rho_{Mg,Mg]')])
latex_labels['MgMg_C'] = OrderedDict([
    ('name',r'C_{Mg,Mg]')])

latex_labels['MgO_A'] = OrderedDict([
    ('name',r'A_{Mg,O}')])
latex_labels['MgO_rho'] = OrderedDict([
    ('name',r'\rho_{Mg,O}')])
latex_labels['MgO_C'] = OrderedDict([
    ('name',r'C_{Mg,O}')])

latex_labels['OO_A'] = OrderedDict([
    ('name',r'A_{OO,OO/}')])
latex_labels['OO_rho'] = OrderedDict([
    ('name',r'\rho_{O,O]')])
latex_labels['OO_C'] = OrderedDict([
    ('name',r'C_{O,O]')])

latex_labels['MgO_NaCl.a0'] = OrderedDict([
    ('name',r'$a_0$'),
    ('units',r'\AA')])
latex_labels['MgO_NaCl.c11'] = OrderedDict([
    ('name',r'$c_{11}$'),
    ('units', r'GPa')])
latex_labels['MgO_NaCl.c12'] = OrderedDict([
    ('name',r'$c_{12}$'),
    ('units', r'GPa')])
latex_labels['MgO_NaCl.c44'] = OrderedDict([
    ('name',r'$c_{44}$'),
    ('units',r'GPa')])
latex_labels['MgO_NaCl.B'] = OrderedDict([
    ('name',r'$B$'),
    ('units',r'GPa')])
latex_labels['MgO_NaCl.G'] = OrderedDict([
    ('name',r'$G$'),
    ('units',r'GPa')])
latex_labels['MgO_NaCl.fr_a'] = OrderedDict([
    ('name',r'$E_{fr,a}$'),
    ('units',r'eV')])
latex_labels['MgO_NaCl.fr_c'] = OrderedDict([
    ('name',r'$E_{fr,c}$'),
    ('units',r'eV')])
latex_labels['MgO_NaCl.sch'] = OrderedDict([
    ('name',r'$E_{sch}$'),
    ('units',r'eV')])
latex_labels['MgO_NaCl.001s'] = OrderedDict([
    ('name',r'$\gamma_{001}$'),
    ('units',r'mJ/m^2')])
latex_labels['MgO_NaCl.optical_ph'] = OrderedDict([
    ('name',r'$\Gamma_{o}$'),
    ('units','cm^-1')])
latex_labels['MgO_NaCl.th_exp'] = OrderedDict([
    ('name',r'$\alpha_{L}$'),
    ('units',r'C^{-1}')])
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
    configuration.qois_validation = qoi_validation_db.qois
    configuration.structures = structure_db
    configuration.potential = potential_formalism
    configuration.sampling_type = sampling
    configuration.sampling_distribution = parameter_distribution
    configuration.sampling_constraints = parameter_constraints
    configuration.reference_potentials = reference_potentials
    configuration.latex_labels = latex_labels
    configuration.write(filename=pyposmat_filename_in)


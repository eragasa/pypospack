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
sampling[0]['n_samples'] = 10000
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'eam'
potential_formalism['symbols'] = ['Ni']
potential_formalism['setfl_filename'] = None
potential_formalism['pair_type'] = 'bornmayer'
potential_formalism['density_type'] = 'eam_dens_exp'
potential_formalism['embedding_type'] = 'eam_embed_eos_rose'

# <---------------- THESE ARE NECESSARY FOR DETERMINING THE SETFL FILE
potential_formalism['N_r'] = 10000
potential_formalism['r_max'] = 10.0
potential_formalism['r_cut'] = 10.0
potential_formalism['N_rho'] = 10000
potential_formalism['rho_max'] = 1000.0
potential_formalism['a0'] = 3.52
potential_formalism['lattice_type'] = 'fcc'

# <---------------- INITIAL PARAMETER DEFINITION
a0 = 3.52
re = 1/(2**0.5)*a0

parameter_distribution = OrderedDict()
parameter_distribution['p_NiNi_phi0'] = [
        'uniform',{
            'a':+0.05,
            'b':+4.00}]
parameter_distribution['p_NiNi_gamma'] = [
        'uniform',{
            'a':2.00,
            'b':7.00}]
parameter_distribution['p_NiNi_r0'] = [
        'equals',re]
parameter_distribution['d_Ni_rho0'] = [
        'uniform',{
            'a':1.0,
            'b':4.0}]
parameter_distribution['d_Ni_beta'] = [
        'uniform',{
            'a':2.0000,
            'b':7.0000}]
parameter_distribution['d_Ni_r0'] = [
        'equals',re]
# parameters - Rose Equation of State
parameter_distribution['e_Ni_latticetype'] = [
        'equals','fcc']
parameter_distribution['e_Ni_ecoh'] = [
        'equals',-4.45]
parameter_distribution['e_Ni_B'] = [
        'equals',188.]
parameter_distribution['e_Ni_a0'] = [
        'equals',a0]
#------------------------------------------------------------------------------
# PARAMETER CONSTRAINTS
#------------------------------------------------------------------------------
parameter_constraints = OrderedDict()
parameter_constraints['p_NiNi_phi0 > 0'] = 'p_NiNi_phi0 > 0.'
parameter_constraints['p_NiNi_gamma > 0'] = 'p_NiNi_gamma > 0.'
parameter_constraints['d_Ni_rho0 > 0'] = 'd_Ni_rho0 > 0.'
parameter_constraints['d_Ni_beta > 0'] = 'd_Ni_beta > 0.'
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['Ni_fcc'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_bcc'] = 'Ni_bcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_sc']  = 'Ni_sc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_hcp'] = 'Ni_hcp_ortho.vasp'
structure_db['structures']['Ni_dia'] = 'Ni_dia_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_fcc_100_unit'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_fcc_110_unit'] = 'Ni_fcc_110_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_fcc_111_unit'] = 'Ni_fcc_111_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_fcc_100_s'] = 'Ni_fcc_100_surf.vasp'
structure_db['structures']['Ni_fcc_110_s'] = 'Ni_fcc_110_surf.vasp'
structure_db['structures']['Ni_fcc_111_s'] = 'Ni_fcc_111_surf.vasp'
structure_db['structures']['Ni_fcc_isf'] = 'Ni_fcc_isf.vasp'
structure_db['structures']['Ni_fcc_esf'] = 'Ni_fcc_esf.vasp'
structure_db['structures']['Ni_fcc_vac'] = 'Ni_fcc_sc_333_vac.vasp'
structure_db['structures']['Ni_fcc_o_int'] = 'Ni_fcc_sc_333_o_int.vasp'
structure_db['structures']['Ni_fcc_i_int'] = 'Ni_fcc_sc_333_t_int.vasp'
#------------------------------------------------------------------------------
# FITTING DATABASE
#------------------------------------------------------------------------------
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='Ni_fcc.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=-4.45)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=3.52)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=261.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=151.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=132.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=188.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=101.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.vac',
        qoi_type='E_formation',
        structures=OrderedDict(
            [
                ('defect','Ni_fcc_vac'),
                ('ideal','Ni_fcc')
            ]
        ),
        target=1.6)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.100s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_100_s'),
                ('ideal','Ni_fcc_100_unit')
            ]
        ),
        target=1.51e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.110s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_110_s'),
                ('ideal','Ni_fcc_110_unit')
            ]
        ),
        target=1.48e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.111s',
        qoi_type='E_surface',
        structures=OrderedDict(
            [
                ('slab','Ni_fcc_111_s'),
                ('ideal','Ni_fcc_111_unit')
            ]
        ),
        target=1.25e-1)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.isf',
        qoi_type='E_stacking_fault',
        structures=OrderedDict([
                ('defect','Ni_fcc_isf'),
                ('ideal','Ni_fcc_111_unit')]),
        target=1.45e-02)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_hcp',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_hcp')]),
        target=0.024)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_bcc',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_bcc')]),
        target=0.092)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_sc',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_sc')]),
        target=0.600)
qoi_db.add_qoi(
        qoi_name='E_Ni_fcc_dia',
        qoi_type='phase_order',
        structures=OrderedDict([
                ('low','Ni_fcc'),
                ('high','Ni_dia')]),
        target=1.27)
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
qoi_constraints['qoi_constraints']=OrderedDict()
qoi_constraints['qoi_constraints']['Ni_fcc.E_coh.abserr'] = ['<',2]
qoi_constraints['qoi_constraints']['Ni_fcc.a0.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.a0']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.c11.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c11']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.c12.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c12']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.c44.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c44']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.B.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c11']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.G.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c12']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.vac'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.c12']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.110s'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.110s']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.100s'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.100s']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.111s'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.111s']['target'])]
qoi_constraints['qoi_constraints']['Ni_fcc.isf.abserr'] = ['<',1.00 * abs(qoi_db.qois['Ni_fcc.isf']['target'])]
qoi_constraints['qoi_constraints']['E_Ni_fcc_bcc'] = ['>',0.]
qoi_constraints['qoi_constraints']['E_Ni_fcc_sc'] = ['>',0.]
qoi_constraints['qoi_constraints']['E_Ni_fcc_hcp'] = ['>',0.]
qoi_constraints['qoi_constraints']['E_Ni_fcc_dia'] = ['>',0.]
qoi_constraints['filter_by__d_zerror'] = OrderedDict()
qoi_constraints['filter_by__d_zerror']['percentile'] = .95
qoi_constraints['select_pareto_only'] = True
#qoi_constraints['filter_by_percentile'] = [80,'pct']
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

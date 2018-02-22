from collections import OrderedDict
from pypospack.qoi import QoiDatabase

a0 = 3.52
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
potential_formalism = OrderedDict()
potential_formalism['potential_type'] = 'eam'
potential_formalism['symbols'] = ['Ni']
potential_formalism['setfl_filename'] = None
potential_formalism['pair_type'] = 'morse'
potential_formalism['density_type'] = 'eam_dens_exp'
potential_formalism['embedding_type'] = 'eam_embed_fs'

# <---------------- THESE ARE NECESSARY FOR DETERMINING THE SETFL FILE
potential_formalism['N_r'] = 10000
potential_formalism['r_max'] = 10.0
potential_formalism['r_cut'] = 10.0
potential_formalism['N_rho'] = 10000
potential_formalism['rho_max'] = 1000.0
potential_formalism['a0'] = a0
potential_formalism['lattice_type'] = 'fcc'

# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
parameter_distribution = OrderedDict()
parameter_distribution['p_NiNi_D0'] = [
        'uniform',{
            'a':0.05,
            'b':0.60}]
parameter_distribution['p_NiNi_a'] = [
        'uniform',{
            'a':0.1,
            'b':4.0}]
parameter_distribution['p_NiNi_r0'] = [
        'uniform',{
            'a':0.1,
            'b':4.0}]
parameter_distribution['d_Ni_rho0'] = [
        'uniform',{
            'a':1.0,
            'b':4.0}]
parameter_distribution['d_Ni_beta'] = [
        'uniform',{
            'a':5.0000,
            'b':7.0000}]
parameter_distribution['d_Ni_r0'] = [
        'uniform',{
            'a':1.00,
            'b':5.00}]
parameter_distribution['e_Ni_F0'] = [
        'uniform',{
            'a':0,
            'b':5.0}]

#------------------------------------------------------------------------------
# PARAMETER CONSTRAINTS
#------------------------------------------------------------------------------
parameter_constraints = OrderedDict()
parameter_constraints['p_NiNi_D0 > 0'] = 'p_NiNi_D0 > 0.'
parameter_constraints['p_NiNi_a > 0'] = 'p_NiNi_a > 0.'
parameter_constraints['p_NiNi_r0 > 0'] = 'p_NiNi_r0 > 0.'
parameter_constraints['d_Ni_rho0 > 0'] = 'd_Ni_rho0 > 0.'
parameter_constraints['d_Ni_beta > 0'] = 'd_Ni_beta > 0.'
parameter_constraints['d_Ni_r0 > 0'] = 'd_Ni_r0 > 0.'
parameter_constraints['e_Ni_F0 > 0'] = 'e_Ni_F0 > 0.'
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['Ni_fcc'] = 'Ni_fcc.vasp'
structure_db['structures']['Ni_fcc_111'] = 'Ni_fcc_111.vasp'
structure_db['structures']['Ni_fcc_110'] = 'Ni_fcc_110.vasp'
structure_db['structures']['Ni_bcc'] = 'Ni_bcc.vasp'
structure_db['structures']['Ni_hcp'] = 'Ni_hcp.vasp'
structure_db['structures']['Ni_sc'] = 'Ni_sc.vasp'
structure_db['structures']['Ni_fcc_100_s'] = 'Ni_fcc_100_s.vasp'
structure_db['structures']['Ni_fcc_110_s'] = 'Ni_fcc_110_s.vasp'
structure_db['structures']['Ni_fcc_111_s'] = 'Ni_fcc_111_s.vasp'
structure_db['structures']['Ni_fcc_usf'] = 'Ni_fcc_usf.vasp'
structure_db['structures']['Ni_fcc_ssf'] = 'Ni_fcc_ssf.vasp'

#------------------------------------------------------------------------------
# FITTING DATABASE 
#------------------------------------------------------------------------------
qoi_db = QoiDatabase()
qoi_db.add_qoi(
        qoi_name='Ni_fcc.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=-5.7771)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=3.508)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=276.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=159.)
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=132.)

qoi_constraints = OrderedDict()
for k,v in qoi_db.qois.items():
    qoi_constraints['{}.err'.format(k)] = 0.20 * abs(qoi_db.qois[k]['target']) 
#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------

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


from collections import OrderedDict
from pypospack.qoi import QoiDatabase

a0 = 3.52
#-----------------------------------------------------------------------------
# DEFINE POTENTIAL FORMALISM
#-----------------------------------------------------------------------------
Ni_eam_potential_formalism = OrderedDict()
Ni_eam_potential_formalism['potential_type'] = 'eam'
Ni_eam_potential_formalism['symbols'] = ['Ni']
Ni_eam_potential_formalism['setfl_filename'] = None
Ni_eam_potential_formalism['pair_type'] = 'morse'
Ni_eam_potential_formalism['density_type'] = 'eam_dens_exp'
Ni_eam_potential_formalism['embedding_type'] = 'eam_embed_universal'

# <---------------- THESE ARE NECESSARY FOR DETERMINING THE SETFL FILE
Ni_eam_potential_formalism['N_r'] = 10000
Ni_eam_potential_formalism['r_max'] = 10.0
Ni_eam_potential_formalism['r_cut'] = 10.0
Ni_eam_potential_formalism['N_rho'] = 10000
Ni_eam_potential_formalism['rho_max'] = 1000.0
Ni_eam_potential_formalism['a0'] = a0
Ni_eam_potential_formalism['lattice_type'] = 'fcc'

r0 = Ni_eam_potential_formalism['a0']/(2**0.5)
# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
Ni_eam_parameter_distribution = OrderedDict()
Ni_eam_parameter_distribution['p_NiNi_D0'] = [
        'uniform',{
            'a':0.050000,
            'b':0.600000}]
Ni_eam_parameter_distribution['p_NiNi_a'] = [
        'uniform',{
            'a':1.000000,
            'b':4.000000}]
Ni_eam_parameter_distribution['p_NiNi_r0'] = [
        'equals',r0]
Ni_eam_parameter_distribution['d_Ni_rho0'] = [
        'uniform',{
            'a':1.000000,
            'b':4.000000}]
Ni_eam_parameter_distribution['d_Ni_beta'] = [
        'uniform',{
            'a':5.0000,
            'b':7.0000}]
Ni_eam_parameter_distribution['d_Ni_r0'] = [
        'equals',r0]
Ni_eam_parameter_distribution['e_Ni_F0'] = [
        'uniform',{
            'a':0.1e-3,
            'b':1.0e-2}]
Ni_eam_parameter_distribution['e_Ni_p'] = [
        'uniform',{
            'a':2.0e-1,
            'b':2.0e+1}]
Ni_eam_parameter_distribution['e_Ni_q'] = [
        'uniform',{
            'a':2.0e-1,
            'b':2.0e+1}]
Ni_eam_parameter_distribution['e_Ni_F1'] = [
        'uniform',{
            'a':-5.,
            'b':+5.}]
Ni_eam_parameter_distribution['e_Ni_rho0'] = [
        'equals',['equilibrium_density',a0,'fcc']
        ]
Ni_eam_parameter_constraints = OrderedDict()
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
Ni_structure_db = OrderedDict()
Ni_structure_db['structure_directory'] = 'structure_db'
Ni_structure_db['structures'] = OrderedDict()
Ni_structure_db['structures']['Ni_fcc'] = 'Ni_fcc.vasp'
Ni_structure_db['structures']['Ni_fcc_111'] = 'Ni_fcc_111.vasp'
Ni_structure_db['structures']['Ni_fcc_110'] = 'Ni_fcc_110.vasp'
Ni_structure_db['structures']['Ni_bcc'] = 'Ni_bcc.vasp'
Ni_structure_db['structures']['Ni_hcp'] = 'Ni_hcp.vasp'
Ni_structure_db['structures']['Ni_sc'] = 'Ni_sc.vasp'
Ni_structure_db['structures']['Ni_fcc_100_s'] = 'Ni_fcc_100_s.vasp'
Ni_structure_db['structures']['Ni_fcc_110_s'] = 'Ni_fcc_110_s.vasp'
Ni_structure_db['structures']['Ni_fcc_111_s'] = 'Ni_fcc_111_s.vasp'
Ni_structure_db['structures']['Ni_fcc_usf'] = 'Ni_fcc_usf.vasp'
Ni_structure_db['structures']['Ni_fcc_ssf'] = 'Ni_fcc_ssf.vasp'

#------------------------------------------------------------------------------
# FITTING DATABASE 
#------------------------------------------------------------------------------
Ni_qoi_db = QoiDatabase()
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.E_coh',
        qoi_type='Ecoh_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=-5.7771)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=3.508)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=276.)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=159.)
Ni_qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=132.)
#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------

# <---------------- SAMPLING CONFIGURATION
Ni_eam_sampling = OrderedDict()
Ni_eam_sampling['n_iterations'] = 30
Ni_eam_sampling['mc_seed'] = 0
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(Ni_eam_sampling['n_iterations']):
    Ni_eam_sampling[i] = OrderedDict()
    Ni_eam_sampling[i]['type'] = 'kde'
    Ni_eam_sampling[i]['n_samples'] = 30
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
Ni_eam_sampling[0]['type'] = 'parametric'


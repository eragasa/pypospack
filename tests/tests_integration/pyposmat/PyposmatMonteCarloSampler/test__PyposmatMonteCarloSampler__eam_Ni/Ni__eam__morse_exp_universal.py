from collections import OrderedDict
from pypospack.qoi import QoiDatabase

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
Ni_eam_potential_formalism['a0'] = 3.52
Ni_eam_potential_formalism['lattice_type'] = 'fcc'

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
        'uniform',{
            'a':2.0,
            'b':3.0}]
Ni_eam_parameter_distribution['d_Ni_rho0'] = [
        'uniform',{
            'a':0.0001,
            'b':20.0000}]
Ni_eam_parameter_distribution['d_Ni_beta'] = [
        'uniform',{
            'a':0.0001,
            'b':20.0000}]
Ni_eam_parameter_distribution['d_Ni_r0'] = [
        'uniform',{
            'a':0.0001,
            'b':10.0000}]
Ni_eam_parameter_distribution['e_Ni_F0'] = [
        'uniform',{
            'a':0.1e-3,
            'b':1.0e-2}]
Ni_eam_parameter_distribution['e_Ni_p'] = [
        'uniform',{
            'a':1.0e-1,
            'b':1.0e+1}]
Ni_eam_parameter_distribution['e_Ni_q'] = [
        'uniform',{
            'a':1.0e-1,
            'b':1.0e+1}]
Ni_eam_parameter_distribution['e_Ni_F1'] = [
        'uniform',{
            'a':-5.,
            'b':+5.}]

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
        target=3.508)
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
Ni_eam_sampling['n_iterations'] = 10
Ni_eam_sampling['mc_seed'] = 0
# <---------------- INITIAL DEFAULT CONFIGURATION
for i in range(Ni_eam_sampling['n_iterations']):
    Ni_eam_sampling[i] = OrderedDict()
    Ni_eam_sampling[i]['type'] = 'kde'
    Ni_eam_sampling[i]['n_samples'] = 10000
# <---------------- OVERRIDE DEFAULT CONFIGURATION, FOR I=0
Ni_eam_sampling[0]['type'] = 'parametric'

#<00--------------- QOI DATABASE
#Si reference database
# LDA-DFT reference data, originally from:
# Pizzagalli et al.  Phil. Mag A (2008) A 83 1191
# Taken from table
# Pizzagalli et al.  J Phys. Condens. Matter (2013) 055801
Si_sw_qoi_db = QoiDatabase()
# <----------------- STRUCTURAL PROPERTIES
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.Ecoh',
    qoi_type='Ecoh_min_all',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=-4.63)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.a0',
    qoi_type='a11_min_all',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=5.43)
# <----------------- ELASTIC PROPERTIES
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c11',
    qoi_type='c11',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=166.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c12',
    qoi_type='c12',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=64.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.c44',
    qoi_type='c44',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=80.0)
Si_sw_qoi_db.add_qoi(
    qoi_name='Si_dia.B',
    qoi_type='bulk_modulus',
    structures=OrderedDict([('ideal','Si_dia')]),
    target=99.0)

# <---------------- QOI PERFORMANCE CONSTRAINTS
# define performance constraints as 20 of the qoi target value
Si_sw_qoi_constraints = OrderedDict()
for qoi_name, qoi_info in Si_sw_qoi_db.qois.items():
    Si_sw_qoi_constraints[qoi_name] = qoi_info['target'] * 0.20

from collections import OrderedDict
from pypospack.qoi import QoiDatabase

#------------------------------------------------------------------------------
# CONFIGURATION SECTION FOR PYPOSMAT PARETO FITTING
#------------------------------------------------------------------------------
# <---------------- SAMPLING CONFIGURATION
sampling = OrderedDict()
sampling['n_iterations'] = 10
sampling['mc_seed'] = None
# <---------------- INITIAL DEFAULT CONFIGURATION
# The default parameter sampler for these runs is sampling for the kernel
# density estimate from the survivors from the previous distribution.
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
potential_formalism['potential_type'] = 'eam'
potential_formalism['symbols'] = ['Ni']
potential_formalism['setfl_filename'] = None
potential_formalism['pair_type'] = 'bornmayer'
potential_formalism['density_type'] = 'eam_dens_exp'
potential_formalism['embedding_type'] = 'eam_embed_fs'

# <---------------- THESE ARE NECESSARY FOR DETERMINING THE SETFL FILE
potential_formalism['N_r'] = 10000
potential_formalism['r_max'] = 10.0
potential_formalism['r_cut'] = 10.0
potential_formalism['N_rho'] = 10000
potential_formalism['rho_max'] = 1000.0
potential_formalism['a0'] = 3.52
potential_formalism['lattice_type'] = 'fcc'

# <---------------- INITIAL PARAMETER DEFINITION
# units need to be in metal units
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
        'uniform',{
            'a':1.00,
            'b':5.00}]
parameter_distribution['d_Ni_rho0'] = [
        'uniform',{
            'a':1.0,
            'b':4.0}]
parameter_distribution['d_Ni_beta'] = [
        'uniform',{
            'a':2.0000,
            'b':7.0000}]
parameter_distribution['d_Ni_r0'] = [
        'uniform',{
            'a':1.00,
            'b':5.00}]
parameter_distribution['e_Ni_F0'] = [
        'uniform',{
            'a':-5.00,
            'b':-0.05}]

#------------------------------------------------------------------------------
# PARAMETER CONSTRAINTS
#------------------------------------------------------------------------------
parameter_constraints = OrderedDict()
parameter_constraints['p_NiNi_phi0 > 0'] = 'p_NiNi_phi0 > 0.'
parameter_constraints['p_NiNi_gamma > 0'] = 'p_NiNi_gamma > 0.'
parameter_constraints['p_NiNi_r0 > 0'] = 'p_NiNi_r0 > 0.'
parameter_constraints['d_Ni_rho0 > 0'] = 'd_Ni_rho0 > 0.'
parameter_constraints['d_Ni_beta > 0'] = 'd_Ni_beta > 0.'
parameter_constraints['d_Ni_r0 > 0'] = 'd_Ni_r0 > 0.'
parameter_constraints['e_Ni_F0 < 0'] = 'e_Ni_F0 < 0.'
#------------------------------------------------------------------------------
# STRUCTURE DATABASE DEFINITION
#------------------------------------------------------------------------------
structure_db = OrderedDict()
structure_db['structure_directory'] = 'structure_db'
structure_db['structures'] = OrderedDict()
structure_db['structures']['Ni_fcc'] = 'Ni_fcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_bcc'] = 'Ni_bcc_100_unit.gga.relaxed.vasp'
structure_db['structures']['Ni_sc']  = 'Ni_sc_100_unit.gga.relaxed.vasp'
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
        target=-4.45 # eV/atom
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.a0',
        qoi_type='a11_min_all',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=3.52 # Angs
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c11',
        qoi_type='c11',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=261. #GPa
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c12',
        qoi_type='c12',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=151. #GPa
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.c44',
        qoi_type='c44',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=132. #GPa
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.B',
        qoi_type='bulk_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=188. #GPa
    )
qoi_db.add_qoi(
        qoi_name='Ni_fcc.G',
        qoi_type='shear_modulus',
        structures=OrderedDict([('ideal','Ni_fcc')]),
        target=101. #GPa
    )
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
qoi_constraints['filter_by_qoi_error']=OrderedDict()
for k,v in qoi_db.qois.items():
    qoi_constraints['filter_by_qoi_error']['{}.err'.format(k)] \
         = 1.00 * abs(qoi_db.qois[k]['target']) 
#qoi_constraints['filter_by_pareto'] = True
qoi_constraints['filter_by_d_metric'] = [80,'pct']
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

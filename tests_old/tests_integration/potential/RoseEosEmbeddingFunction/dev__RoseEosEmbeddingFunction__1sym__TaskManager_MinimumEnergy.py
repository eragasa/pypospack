import os
from pypospack.task import TaskManager
from collections import OrderedDict
pyposmat_root_dir = [v.strip() for v in os.environ['PYTHONPATH'].split(':') if v.endswith('pypospack')][0]
_base_directory = "./"

pot_definition = OrderedDict()
pot_definition['potential_type'] = 'eam'
pot_definition['symbols'] = ['Ni']
pot_definition['setfl_filename'] = None
pot_definition['pair_type'] = 'bornmayer'
pot_definition['density_type'] = 'eam_dens_exp'
pot_definition['embedding_type'] = 'eam_embed_eos_rose'
pot_definition['N_r'] = 10000
pot_definition['r_max'] = 10.0
pot_definition['r_cut'] = 10.0
pot_definition['N_rho'] = 10000
pot_definition['rho_max'] = 1000.0
pot_definition['symbols'] = ['Ni']
pot_definition['a0'] = 3.52
pot_definition['lattice_type'] = 'fcc'

_tasks = OrderedDict([
    ('Ni_fcc.lmps_min_all', OrderedDict([
        ('task_type', 'lmps_min_all'), 
        ('task_structure', 'Ni_fcc')])), 
#    ('Ni_fcc.lmps_elastic', OrderedDict([
#        ('task_type', 'lmps_elastic'), 
#        ('task_structure', 'Ni_fcc')])), 
#    ('Ni_fcc_vac.lmps_min_pos', OrderedDict([
#        ('task_type', 'lmps_min_pos'), 
#        ('task_structure', 'Ni_fcc_vac'), 
#        ('bulk_structure', 'Ni_fcc')])), 
#    ('Ni_fcc_100_unit.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_fcc_100_unit')])), 
#    ('Ni_fcc_100_s.lmps_min_pos', OrderedDict([
#        ('task_type', 'lmps_min_pos'), 
#        ('task_structure', 'Ni_fcc_100_s'), 
#        ('bulk_structure', 'Ni_fcc_100_unit')])), 
#    ('Ni_fcc_110_unit.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_fcc_110_unit')])), 
#    ('Ni_fcc_110_s.lmps_min_pos', OrderedDict([
#        ('task_type', 'lmps_min_pos'), 
#        ('task_structure', 'Ni_fcc_110_s'), 
#        ('bulk_structure', 'Ni_fcc_110_unit')])), 
#    ('Ni_fcc_111_unit.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_fcc_111_unit')])), 
#    ('Ni_fcc_111_s.lmps_min_pos', OrderedDict([
#        ('task_type', 'lmps_min_pos'), 
#        ('task_structure', 'Ni_fcc_111_s'), 
#        ('bulk_structure', 'Ni_fcc_111_unit')])), 
#    ('Ni_fcc_isf.lmps_min_sf', OrderedDict([
#        ('task_type', 'lmps_min_sf'), 
#        ('task_structure', 'Ni_fcc_isf'), 
#        ('bulk_structure', 'Ni_fcc_111_unit')])), 
#    ('Ni_hcp.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_hcp')])), 
#    ('Ni_bcc.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_bcc')])), 
#    ('Ni_sc.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_sc')])), 
#    ('Ni_dia.lmps_min_all', OrderedDict([
#        ('task_type', 'lmps_min_all'), 
#        ('task_structure', 'Ni_dia')]))
    ])


_structures = OrderedDict([
    ('structure_directory', os.path.join(pyposmat_root_dir,'data/Ni_structure_db')),
    ('structures', OrderedDict([
        ('Ni_fcc', 'Ni_fcc_100_unit.gga.relaxed.vasp'), 
        ('Ni_bcc', 'Ni_bcc_100_unit.gga.relaxed.vasp'), 
        ('Ni_sc', 'Ni_sc_100_unit.gga.relaxed.vasp'), 
        ('Ni_hcp', 'Ni_hcp_ortho.vasp'), 
        ('Ni_dia', 'Ni_dia_100_unit.gga.relaxed.vasp'), 
        ('Ni_fcc_100_unit', 'Ni_fcc_100_unit.gga.relaxed.vasp'), 
        ('Ni_fcc_110_unit', 'Ni_fcc_110_unit.gga.relaxed.vasp'), 
        ('Ni_fcc_111_unit', 'Ni_fcc_111_unit.gga.relaxed.vasp'), 
        ('Ni_fcc_100_s', 'Ni_fcc_100_surf.vasp'), 
        ('Ni_fcc_110_s', 'Ni_fcc_110_surf.vasp'), 
        ('Ni_fcc_111_s', 'Ni_fcc_111_surf.vasp'), 
        ('Ni_fcc_isf', 'Ni_fcc_isf.vasp'), 
        ('Ni_fcc_esf', 'Ni_fcc_esf.vasp'), 
        ('Ni_fcc_vac', 'Ni_fcc_sc_333_vac.vasp'), 
        ('Ni_fcc_o_int', 'Ni_fcc_sc_333_o_int.vasp'), 
        ('Ni_fcc_i_int', 'Ni_fcc_sc_333_t_int.vasp')]))
    ])


a0=3.52
r0=1/(2**0.5)*a0

parameters = OrderedDict()
parameters['p_NiNi_phi0'] = 1.0
parameters['p_NiNi_gamma'] = 2.0
parameters['p_NiNi_r0'] = r0
parameters['d_Ni_rho0'] = 1.0
parameters['d_Ni_beta'] = 4.0
parameters['d_Ni_r0'] = r0
parameters['e_Ni_latticetype'] = 'fcc'
parameters['e_Ni_ecoh'] = -4.45
parameters['e_Ni_B']= 188.
parameters['e_Ni_a0'] = a0

task_manager = TaskManager(base_directory=_base_directory)
task_manager.configure(
        tasks=_tasks,
        structures=_structures)
results = task_manager.evaluate_tasks(
        parameters=parameters,
        potential=pot_definition)
print(results)

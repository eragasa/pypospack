import os,shutil
from collections import OrderedDict
from pypospack.task.lammps import LammpsStructuralMinimization

Ni_eam_1_potential_definition = OrderedDict()
Ni_eam_1_potential_definition['potential_type'] = 'eam'
Ni_eam_1_potential_definition['symbols'] = ['Ni']
Ni_eam_1_potential_definition['setfl_filename'] = None
Ni_eam_1_potential_definition['pair_type'] = 'morse'
Ni_eam_1_potential_definition['density_type'] = 'eam_dens_exp'
Ni_eam_1_potential_definition['embedding_type'] = 'eam_embed_universal'
Ni_eam_1_potential_definition['N_r'] = 10000
Ni_eam_1_potential_definition['r_max'] = 10.0
Ni_eam_1_potential_definition['r_cut'] = 10.0
Ni_eam_1_potential_definition['N_rho'] = 10000
Ni_eam_1_potential_definition['rho_max'] = 1000.0
Ni_eam_1_potential_definition['symbols'] = ['Ni']
Ni_eam_1_potential_definition['a0'] = 3.52
Ni_eam_1_potential_definition['lattice_type'] = 'fcc'

Ni_eam_1_parameters = OrderedDict()
Ni_eam_1_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_1_parameters['p_NiNi_a'] = 3.429506
Ni_eam_1_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_1_parameters['d_Ni_rho0'] = 10.0
Ni_eam_1_parameters['d_Ni_beta'] = 5.0
Ni_eam_1_parameters['d_Ni_r0'] = 2.0
Ni_eam_1_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_1_parameters['e_Ni_p'] = 8.96274624
Ni_eam_1_parameters['e_Ni_q'] = 8.95940869
Ni_eam_1_parameters['e_Ni_F1'] = -1.09

Ni_structures = OrderedDict()
Ni_structures['structure_directory'] = 'structure_db'
Ni_structures['structures'] = OrderedDict()
Ni_structures['structures']['Ni_fcc'] = 'Ni_fcc_unit_001.vasp'

Ni_task_configuration = OrderedDict()
Ni_task_configuration['task'] = OrderedDict()
Ni_task_configuration['task']['task_name'] = 'Ni_fcc_unit.lmps_min_all'
Ni_task_configuration['task']['task_directory'] = 'Ni_fcc_unit.lmps_min_all'
Ni_task_configuration['task']['task_type'] = 'lmps_min_all'
Ni_task_configuration['potential'] = Ni_eam_1_potential_definition
Ni_task_configuration['parameters'] = Ni_eam_1_parameters
Ni_task_configuration['structure'] = Ni_structures

configuration = Ni_task_configuration

task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = os.path.join(
        configuration['structure']['structure_directory'],
        configuration['structure']['structures']['Ni_fcc']
    )
restart=False
fullauto=False

lammps_task = LammpsStructuralMinimization(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)
lammps_task.on_running(configuration)
while lammps_task.status is not 'POST':
    lammps_task.update_status()
lammps_task.on_post(configuration)

print(lammps_task.results)

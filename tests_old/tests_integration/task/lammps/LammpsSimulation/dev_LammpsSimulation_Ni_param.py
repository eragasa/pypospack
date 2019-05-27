import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

Ni_task_information = OrderedDict()
Ni_task_information['task_name'] = 'Ni_fcc_unit.E_min_all'
Ni_task_information['task_directory'] = 'Ni_fcc_unit.E_min_all'
Ni_task_information['task_type'] = 'lammps_relax_none'

Ni_eam_potential_definition = OrderedDict()
Ni_eam_potential_definition['potential_type'] = 'eam'
Ni_eam_potential_definition['setfl_filename']=None
Ni_eam_potential_definition['pair_type']='morse'
Ni_eam_potential_definition['density_type']='eam_dens_exp'
Ni_eam_potential_definition['embedding_type']='eam_embed_universal'
Ni_eam_potential_definition['N_r'] = 2000
Ni_eam_potential_definition['r_max'] = 10.0
Ni_eam_potential_definition['r_cut'] = 10.0
Ni_eam_potential_definition['N_rho'] = 2000
Ni_eam_potential_definition['rho_max'] = 10.0
Ni_eam_potential_definition['symbols'] = ['Ni']

Ni_eam_parameters = OrderedDict()
Ni_eam_parameters['p_NiNi_D0'] = 0.001114
Ni_eam_parameters['p_NiNi_a'] = 3.429506
Ni_eam_parameters['p_NiNi_r0'] = 2.6813
Ni_eam_parameters['d_Ni_rho0'] = 10.0
Ni_eam_parameters['d_Ni_beta'] = 5.0
Ni_eam_parameters['d_Ni_r0'] = 2.0
Ni_eam_parameters['e_Ni_F0'] = 4.10341782e-3
Ni_eam_parameters['e_Ni_p'] = 8.96274624
Ni_eam_parameters['e_Ni_q'] = 8.95940869
Ni_eam_parameters['e_Ni_F1'] = -3.09

Ni_structure_definition = OrderedDict()
Ni_structure_definition['name'] = 'Ni_fcc_unit'
Ni_structure_definition['filename'] = os.path.join(
        'test_LammpsSimulation',
        'Ni_fcc_unit.vasp')

Ni_task_configuration= OrderedDict()
Ni_task_configuration['task'] = OrderedDict()
Ni_task_configuration['task']['task_name'] = 'Ni_fcc_unit.E_min_all'
Ni_task_configuration['task']['task_directory'] = 'Ni_fcc_unit.E_min_all'
Ni_task_configuration['task']['task_type'] = 'single_point'
Ni_task_configuration['potential'] = Ni_eam_potential_definition
Ni_task_configuration['parameters'] = Ni_eam_parameters
Ni_task_configuration['structure'] = Ni_structure_definition

configuration=Ni_task_configuration

task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = configuration['structure']['filename']
restart=False
fullauto=False

from pypospack.task.lammps import LammpsSimulation
lammps_task = LammpsSimulation(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)
lammps_task.on_running(configuration)
lammps_task.on_post(configuration)

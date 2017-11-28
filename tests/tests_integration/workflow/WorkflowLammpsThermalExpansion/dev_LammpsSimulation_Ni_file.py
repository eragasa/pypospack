import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

Ni_eam_potential_definition = OrderedDict()
Ni_eam_potential_definition['potential_type'] = 'eam'
Ni_eam_potential_definition['symbols'] = ['Ni']
Ni_eam_potential_definition['setfl_filename'] = os.path.join(
        'test_LammpsSimulation',
        'Ni_Mendelev_2010.eam.fs')

print('setfl_filename:{}'.format(Ni_eam_potential_definition['setfl_filename']))
assert os.path.isfile(Ni_eam_potential_definition['setfl_filename'])

Ni_eam_parameters = None

Ni_structure_definition = OrderedDict()
Ni_structure_definition['name'] = 'Ni_fcc_unit'
Ni_structure_definition['filename'] = os.path.join(
        'test_LammpsSimulation',
        'Ni_fcc_unit.vasp')

print('poscar_filename={}'.format(Ni_structure_definition['filename']))
assert os.path.isfile(Ni_structure_definition['filename'])

Ni_task_configuration= OrderedDict()
Ni_task_configuration['task'] = OrderedDict()
Ni_task_configuration['task']['task_name'] = 'Ni_fcc_unit.E_min_all'
Ni_task_configuration['task']['task_directory'] = 'Ni_fcc_.E_min_all'
Ni_task_configuration['task_type'] = 'min_none'
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

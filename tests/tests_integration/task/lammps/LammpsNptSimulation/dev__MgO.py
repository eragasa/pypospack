import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

MgO_buck_potential_definition = OrderedDict()
MgO_buck_potential_definition['potential_type'] = 'buckingham'
MgO_buck_potential_definition['symbols'] = ['Mg','O']

MgO_LC_parameters = OrderedDict()
MgO_LC_parameters['chrg_Mg'] = +2.0
MgO_LC_parameters['chrg_O']  = -2.0
MgO_LC_parameters['MgMg_A']   = 0.0 
MgO_LC_parameters['MgMg_rho'] = 0.5
MgO_LC_parameters['MgMg_C']   = 0.0
MgO_LC_parameters['MgO_A']    = 821.6
MgO_LC_parameters['MgO_rho']  = 0.3242
MgO_LC_parameters['MgO_C']    = 0.0
MgO_LC_parameters['OO_A']     = 2274.00 
MgO_LC_parameters['OO_rho']   = 0.1490
MgO_LC_parameters['OO_C']     = 27.88

MgO_structure_definition = OrderedDict()
MgO_structure_definition['name'] = 'MgO_NaCl_unit'
MgO_structure_definition['filename'] = os.path.join(
        '../../../../structure_db/MgO_structure_db/',
        'MgO_NaCl_unit.gga.relax.vasp')

MgO_LC_configuration = OrderedDict()
MgO_LC_configuration['task'] = OrderedDict()
MgO_LC_configuration['task']['task_name'] = 'MgO_NaCl.themrmal_expansion'
MgO_LC_configuration['task']['task_directory'] = 'MgO_NaCl.thermal_expansion'
MgO_LC_configuration['task_type'] = 'min_none'
MgO_LC_configuration['potential'] = MgO_buck_potential_definition
MgO_LC_configuration['parameters'] = MgO_LC_parameters
MgO_LC_configuration['structure'] = MgO_structure_definition

configuration=MgO_LC_configuration

task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = configuration['structure']['filename']
restart=False
fullauto=False

from pypospack.task.lammps import LammpsNptSimulation
lammps_task = LammpsNptSimulation(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename,
        temperature = 100,
        time_total = 10,
        time_step = 0.001,
        supercell = [5,5,5],
        pressure = 0)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)

while lammps_task.status is 'RUNNING':
    lammps_task.on_running(configuration)

lammps_task.on_post(configuration)

print(lammps_task.results)

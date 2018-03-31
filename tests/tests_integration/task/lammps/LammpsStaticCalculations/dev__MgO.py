import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

potential_definition = OrderedDict()
potential_definition['potential_type'] = 'buckingham'
potential_definition['symbols'] = ['Mg','O']

parameters = OrderedDict()
parameters['chrg_Mg'] = +2.0
parameters['chrg_O']  = -2.0
parameters['MgMg_A']   = 0.0 
parameters['MgMg_rho'] = 0.5
parameters['MgMg_C']   = 0.0
parameters['MgO_A']    = 821.6
parameters['MgO_rho']  = 0.3242
parameters['MgO_C']    = 0.0
parameters['OO_A']     = 2274.00 
parameters['OO_rho']   = 0.1490
parameters['OO_C']     = 27.88

strucutre_db = OrderedDict()
strucutre_db['name'] = 'MgO_NaCl_unit'
strucutre_db['filename'] = os.path.join(
        'structure_db',
        'MgO_NaCl_unit.gga.relax.vasp')

configuration = OrderedDict()
configuration['task'] = OrderedDict()
configuration['task']['task_name'] = 'MgO_NaCl.E_sp'
configuration['task']['task_directory'] = 'MgO_NaCl.E_sp'
configuration['task_type'] = 'min_none'
configuration['potential'] = potential_definition
configuration['parameters'] = parameters
configuration['structure'] = strucutre_db


task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = configuration['structure']['filename']
restart=False
fullauto=False

from pypospack.task.lammps import LammpsStaticCalculations
lammps_task = LammpsStaticCalculations(
        task_name = task_name,
        task_directory = task_directory,
        structure_filename = structure_filename)
lammps_task.on_init(configuration)
lammps_task.on_config(configuration)
lammps_task.on_ready(configuration)
lammps_task.on_running(configuration)
while lammps_task.status != 'POST':
    lammps_task.update_status()
lammps_task.on_post(configuration)

print(lammps_task.results)

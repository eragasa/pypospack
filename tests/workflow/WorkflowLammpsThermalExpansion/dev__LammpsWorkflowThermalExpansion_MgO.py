import os,shutil,time
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential
from pypospack.task.lammps import LammpsNptSimulation

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
        'test_WorkflowLammpsThermalExpansion',
        'MgO_NaCl_unit.gga.relax.vasp')
MgO_structure_definition['supercell'] = [10,10,10]

MgO_LC_configuration = OrderedDict()
MgO_LC_configuration['task'] = OrderedDict()
MgO_LC_configuration['task']['task_name'] = 'MgO_NaCl_unit.npt.T1000'
MgO_LC_configuration['task']['task_directory'] = 'MgO_NaCl_unit.npt.T1000'
MgO_LC_configuration['task']['dt'] = 0.001 # in picoseconds
MgO_LC_configuration['task']['temperature'] = 1000 # in Kelvin
MgO_LC_configuration['task_type'] = 'lmps_npt'
MgO_LC_configuration['potential'] = MgO_buck_potential_definition
MgO_LC_configuration['parameters'] = MgO_LC_parameters
MgO_LC_configuration['structure'] = MgO_structure_definition

configuration=MgO_LC_configuration

temp_0 = 100
temp_f = 2000
for temperature in range(temp_0,temp_f+1,100):
    configuration['task']['task_name'] = 'MgO_NaCl_unit.npt.T{}'.format(
            temperature)
    configuration['task']['task_directory'] = 'MgO_NaCl_unit.npt.T{}'.format(
            temperature)
    configuration['task']['temperature'] = temperature
    configuration['task']['task_type'] = 'lmps_npt'

    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)
    
    lammps_task = LammpsNptSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    lammps_task.on_ready(configuration)

    _time_sleep = 0.1
    time.sleep(_time_sleep)
    print("t={}".format(temperature))
    while (lammps_task.process.poll() is not None):
        time.sleep(_time_sleep)
    lammps_task = None


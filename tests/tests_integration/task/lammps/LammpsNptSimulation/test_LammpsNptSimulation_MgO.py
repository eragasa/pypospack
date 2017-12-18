import pytest
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
        'test_LammpsNptSimulation',
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

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsNptSimulation

def test____init___():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    #<--- code being tested
    from pypospack.task.lammps import LammpsNptSimulation
    lammps_task = LammpsNptSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)

    #<--- check directory structure
    assert os.path.isdir(
            os.path.abspath(lammps_task.task_directory))
    assert os.listdir(lammps_task.task_directory) == []
    #<--- check attributes
    assert lammps_task.task_name == task_name
    assert os.path.abspath(lammps_task.task_directory)\
            == os.path.abspath(task_directory)
    assert lammps_task.task_type == 'single_point'
    assert lammps_task.lammps_bin == os.environ['LAMMPS_BIN']
    assert lammps_task.lammps_input_filename == 'lammps.in'
    assert lammps_task.lammps_output_filename == 'lammps.out'
    assert lammps_task.lammps_structure_filename == 'lammps.structure'
    assert lammps_task.lammps_setfl_filename is None
    assert lammps_task.potential is None
    assert lammps_task.structure_filename == structure_filename
    assert isinstance(lammps_task.structure,crystal.SimulationCell)

    assert lammps_task.status == 'INIT'

def test__on_init():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    assert 'potential' in configuration
    assert 'parameters' in configuration
    #<--- code setup
    from pypospack.task.lammps import LammpsNptSimulation
    lammps_task = LammpsNptSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'INIT'
    assert lammps_task.potential is None
    #<--- code being testing
    lammps_task.on_init(configuration)

    #<--- check directory structure
    assert os.path.isdir(
            os.path.abspath(lammps_task.task_directory))
    #<--- check attributes
    assert lammps_task.task_name == task_name
    assert os.path.abspath(lammps_task.task_directory)\
            == os.path.abspath(task_directory)
    assert lammps_task.task_type == 'single_point'
    assert lammps_task.lammps_bin == os.environ['LAMMPS_BIN']
    assert lammps_task.lammps_input_filename == 'lammps.in'
    assert lammps_task.lammps_output_filename == 'lammps.out'
    assert lammps_task.lammps_structure_filename == 'lammps.structure'
    assert lammps_task.lammps_setfl_filename is None
    assert isinstance(lammps_task.potential,potential.Potential)
    assert lammps_task.structure_filename == structure_filename
    assert isinstance(lammps_task.structure,crystal.SimulationCell)
    assert lammps_task.process is None

def test__on_ready():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    assert 'potential' in configuration
    assert 'parameters' in configuration
    #<--- code setup
    from pypospack.task.lammps import LammpsNptSimulation
    lammps_task = LammpsNptSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'READY'
    #<--- code being testing
    lammps_task.on_ready(configuration)

    assert type(lammps_task.conditions_INIT) == OrderedDict
    assert type(lammps_task.conditions_CONFIG) == OrderedDict
    assert type(lammps_task.conditions_READY) == OrderedDict
    assert type(lammps_task.conditions_RUNNING) == OrderedDict
    assert type(lammps_task.conditions_POST) == OrderedDict
    assert type(lammps_task.conditions_FINISHED) == OrderedDict

    assert lammps_task.conditions_INIT['task_directory_created']
    assert all([v for k,v in lammps_task.conditions_INIT.items()]) == True
    assert lammps_task.conditions_CONFIG['potential_initialized'] == True
    assert lammps_task.conditions_CONFIG['parameters_processed'] == True
    assert all([v for k,v in lammps_task.conditions_CONFIG.items()]) == True
    assert all([v for k,v in lammps_task.conditions_READY.items()]) == True
    assert lammps_task.conditions_RUNNING['process_initialized']== True
    assert all([v for k,v in lammps_task.conditions_RUNNING.items()]) == True
    assert lammps_task.conditions_POST['process_finished'] == False
    assert all([v for k,v in lammps_task.conditions_POST.items()]) == False
    assert all([v for k,v in lammps_task.conditions_FINISHED.items()]) == False


import pytest
import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

potential_definition = OrderedDict()
potential_definition['potential_type'] = 'buckingham'
potential_definition['symbols'] = ['Mg','O']

# parameters from Lewis and Catlow
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

structure_db = OrderedDict()
structure_db['name'] = 'MgO_NaCl_unit'
structure_db['filename'] = os.path.join(
        'structure_db',
        'MgO_NaCl_unit.gga.relax.vasp')

configuration = OrderedDict()
configuration['task'] = OrderedDict()
configuration['task']['task_name'] = 'MgO_NaCl.E_min_none'
configuration['task']['task_directory'] = 'MgO_NaCl.E_min_none'
configuration['task']['task_type'] = 'lmps_min_none'
configuration['potential'] = potential_definition
configuration['parameters'] = parameters
configuration['structure'] = structure_db

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsStaticCalculations

def test____init___():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    #<--- code being tested
    from pypospack.task.lammps import LammpsStaticCalculations
    lammps_task = LammpsStaticCalculations(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    #<--- expected behavior
    task_type = configuration['task']['task_type']
    #<--- check directory structure
    assert os.path.isdir(
            os.path.abspath(lammps_task.task_directory))
    assert os.listdir(lammps_task.task_directory) == []
    #<--- check attributes
    assert lammps_task.task_name == task_name
    assert os.path.abspath(lammps_task.task_directory)\
            == os.path.abspath(task_directory)
    assert lammps_task.task_type == task_type
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
    from pypospack.task.lammps import LammpsStaticCalculations
    lammps_task = LammpsStaticCalculations(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'INIT'
    assert lammps_task.potential is None
    #<--- code being testing
    lammps_task.on_init(configuration)
    #<--- expected behavior
    task_type = configuration['task']['task_type']
    #<--- check directory structure
    assert os.path.isdir(
            os.path.abspath(lammps_task.task_directory))
    #<--- check attributes
    assert lammps_task.task_name == task_name
    assert os.path.abspath(lammps_task.task_directory)\
            == os.path.abspath(task_directory)
    assert lammps_task.task_type == task_type
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
    from pypospack.task.lammps import LammpsStaticCalculations
    lammps_task = LammpsStaticCalculations(
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

    #if len(lammps_task.conditions_READY) == 0:
    #    assert lammps_task.status == 'READY'
    #else:
    #    assert lammps_task.status == 'CONFIG'


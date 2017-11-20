import pytest
import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

Ni_eam_potential_definition = OrderedDict()
Ni_eam_potential_definition['potential_type'] = 'eam'
Ni_eam_potential_definition['symbols'] = ['Ni']
Ni_eam_potential_definition['setfl_filename'] = os.path.join(
        'test_LammpsSinglePointCalculation',
        'Ni_Mendelev_2010.eam.fs')

print(Ni_eam_potential_definition['setfl_filename'])
assert os.path.isfile(Ni_eam_potential_definition['setfl_filename'])

Ni_eam_parameters = None

Ni_structure_definition = OrderedDict()
Ni_structure_definition['name'] = 'Ni_fcc_unit'
Ni_structure_definition['filename'] = os.path.join(
        'test_LammpsSinglePointCalculation',
        'Ni_fcc_unit.vasp')

print(Ni_structure_definition['filename'])
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

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsSinglePointCalculation

def test____init___():
    symbols = configuration['potential']['symbols']
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    eam_setfl_filename = configuration['potential']['setfl_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    assert type(eam_setfl_filename) is str
    assert type(structure_filename) is str
    poscar = vasp.Poscar()
    poscar.read(structure_filename)
    assert os.path.isfile(eam_setfl_filename)
    assert os.path.isfile(structure_filename)
    #<--- code being tested
    from pypospack.task.lammps import LammpsSinglePointCalculation
    lammps_task = LammpsSinglePointCalculation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)

    #<--- expected behavior
    lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))
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
    assert lammps_task.lammps_setfl_filename == None
    assert lammps_task.potential is None
    assert lammps_task.structure_filename == structure_filename
    assert isinstance(lammps_task.structure,crystal.SimulationCell)

    assert lammps_task.status == 'INIT'

def test__on_init():
    symbols = configuration['potential']['symbols']
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
    from pypospack.task.lammps import LammpsSinglePointCalculation
    lammps_task = LammpsSinglePointCalculation(
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
    lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))
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
    symbols = configuration['potential']['symbols']
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
    from pypospack.task.lammps import LammpsSinglePointCalculation
    lammps_task = LammpsSinglePointCalculation(
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
    #<--- expected behavior
    lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))

    #<--- check to see if the configuration dictionaries were setup correctly
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
    assert lammps_task.conditions_RUNNING['process_initialized']== False
    assert all([v for k,v in lammps_task.conditions_RUNNING.items()]) == False
    assert lammps_task.conditions_POST['process_finished'] == False
    assert all([v for k,v in lammps_task.conditions_POST.items()]) == False
    assert all([v for k,v in lammps_task.conditions_FINISHED.items()]) == False

    if len(lammps_task.conditions_READY) == 0:
        assert lammps_task.status == 'READY'
    else:
        assert lammps_task.status == 'CONFIG'


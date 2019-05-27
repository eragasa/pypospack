import pytest
import os,shutil,time
from collections import OrderedDict
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
from pypospack.potential import Potential

from MgO_buck import MgO_LewisCatlow
from MgO_buck import MgO_structures

MgO_LC_configuration = OrderedDict()
MgO_LC_configuration['task'] = OrderedDict()
MgO_LC_configuration['task']['task_name'] = 'MgO_NaCl.lmps_min_all'
MgO_LC_configuration['task']['task_directory'] = 'MgO_NaCl.lmps_min_all'
MgO_LC_configuration['task']['structure_filename'] = os.path.join(
        MgO_structures['structure_db_dir'],
        MgO_structures['MgO_NaCl_unit']['filename'])

MgO_LC_configuration['potential'] = MgO_LewisCatlow['potential']
MgO_LC_configuration['parameters'] = MgO_LewisCatlow['parameters']

MgO_LC_configuration['expected'] = OrderedDict()
MgO_LC_configuration['expected']['task_type'] = 'lmps_min_all'

configuration=MgO_LC_configuration

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsStructuralMinimization

def test____init___():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['task']['structure_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    #<--- code being tested
    from pypospack.task.lammps import LammpsStructuralMinimization
    lammps_task = LammpsStructuralMinimization(
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
    structure_filename = configuration['task']['structure_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    assert 'potential' in configuration
    assert 'parameters' in configuration
    #<--- code setup
    from pypospack.task.lammps import LammpsStructuralMinimization
    lammps_task = LammpsStructuralMinimization(
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
    assert isinstance(lammps_task.potential,Potential)
    assert lammps_task.structure_filename == structure_filename
    assert isinstance(lammps_task.structure,crystal.SimulationCell)
    assert lammps_task.process is None

def test__on_ready():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['task']['structure_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    assert 'potential' in configuration
    assert 'parameters' in configuration
    #<--- code setup
    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
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
    #assert lammps_task.conditions_RUNNING['process_initialized']== False
    #assert all([v for k,v in lammps_task.conditions_RUNNING.items()]) == False
    assert lammps_task.conditions_POST['process_finished'] == False
    assert all([v for k,v in lammps_task.conditions_POST.items()]) == False
    assert all([v for k,v in lammps_task.conditions_FINISHED.items()]) == False

    #if len(lammps_task.conditions_READY) == 0:
    #    assert lammps_task.status == 'READY'
    #else:
    #    assert lammps_task.status == 'CONFIG'

def test__on_post():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['task']['structure_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    #<--- code setup
    from pypospack.task.lammps import LammpsStructuralMinimization
    lammps_task = LammpsStructuralMinimization(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    lammps_task.on_ready(configuration)
    while lammps_task.status != 'POST':
        lammps_task.update_status()
        time.sleep(0.1)
    lammps_task.on_post(configuration)
    assert isinstance(lammps_task.results,OrderedDict)


import pytest
import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

# An example script for running a LAMMPS simulation for the Stillinger-Weber 
# parameterization
# References:
#     http://lammps.sandia.gov/doc/pair_sw.html

Si_sw_potential_definition = OrderedDict()
Si_sw_potential_definition['potential_type'] = 'stillingerweber'
Si_sw_potential_definition['symbols'] = ['Si']

# Potential parameters from Pizzagalli, J Phys Condens Matter 25 (2013) 055801
Si_sw_pizzagalli = OrderedDict()
Si_sw_pizzagalli['SiSiSi_epsilon']=1.04190
Si_sw_pizzagalli['SiSiSi_sigma']=2.128117
Si_sw_pizzagalli['SiSiSi_a']=1.80
Si_sw_pizzagalli['SiSiSi_lambda']=31.0
Si_sw_pizzagalli['SiSiSi_gamma']=1.10
Si_sw_pizzagalli['SiSiSi_costheta0']=-1/3
Si_sw_pizzagalli['SiSiSi_A']=19.0
Si_sw_pizzagalli['SiSiSi_B']=0.65
Si_sw_pizzagalli['SiSiSi_p']=3.5
Si_sw_pizzagalli['SiSiSi_q']=0.5
Si_sw_pizzagalli['SiSiSi_tol']=0

src_structure_dir = 'test_LammpsElasticCalculation'
Si_structure_definition = OrderedDict()
Si_structure_definition['name'] = 'Si_sc'
Si_structure_definition['filename'] = os.path.join(
        src_structure_dir,
        'Si_sc_unit.vasp')

Si_configuration = OrderedDict()
Si_configuration['task'] = OrderedDict()
Si_configuration['task']['task_name'] = 'Si_sc_unit.lmps_min_all'
Si_configuration['task']['task_directory'] = 'Si_sc_unit.lmps_min_all'
Si_configuration['task']['task_type'] = 'lmps_min_all'
Si_configuration['potential'] = Si_sw_potential_definition
Si_configuration['parameters'] = Si_sw_pizzagalli
Si_configuration['structure'] = Si_structure_definition

configuration=Si_configuration

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsElasticCalculation

def test____init___():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    #<--- code being tested
    from pypospack.task.lammps import LammpsElasticCalculation
    lammps_task = LammpsElasticCalculation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)

    lammps_bin = os.environ['LAMMPS_BIN']
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
    assert lammps_task.lammps_bin == lammps_bin
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
    from pypospack.task.lammps import LammpsElasticCalculation
    lammps_task = LammpsElasticCalculation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'INIT'
    assert lammps_task.potential is None
    #<--- code being testing
    lammps_task.on_init(configuration)
    #<--- expected results
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

    if len(lammps_task.conditions_READY) == 0:
        assert lammps_task.status == 'READY'
    else:
        assert lammps_task.status == 'CONFIG'

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
    from pypospack.task.lammps import LammpsElasticCalculation
    lammps_task = LammpsElasticCalculation(
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

def test__on_post():
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
    from pypospack.task.lammps import LammpsElasticCalculation
    lammps_task = LammpsElasticCalculation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.on_init(configuration)
    lammps_task.on_config(configuration)
    lammps_task.on_ready(configuration)
    #<--- code being testing
    while lammps_task.status != 'POST':
        lammps_task.update_status()
    lammps_task.on_post(configuration)
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'READY'

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


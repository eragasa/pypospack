import pytest
import os,shutil
from collections import OrderedDict

import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

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

assert configuration['potential']['potential_type']

task_name = configuration['task']['task_name']
task_directory = configuration['task']['task_directory']
structure_filename = configuration['structure']['filename']
restart=False
fullauto=False

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsSimulation

def validate_configuration(configuration):
    assert not os.path.exists(task_directory)

    eam_setfl_filename = configuration['potential']['setfl_filename']
    if type(eam_setfl_filename) is str:
        assert os.path.isfile(eam_setfl_filename)
    else:
        assert eam_setfl_filename is None

    structure_filename = configuration['structure']['filename']
    assert type(structure_filename) is str
    assert os.path.isfile(structure_filename)
    poscar = vasp.Poscar()
    poscar.read(structure_filename)

def test____init___():
    symbols = configuration['potential']['symbols']
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    eam_setfl_filename = configuration['potential']['setfl_filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    validate_configuration(configuration) 
    #<--- code being tested
    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)

    #<--- expected behavior
    lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))
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
    assert lammps_task.lammps_setfl_filename == None
    assert lammps_task.potential is None
    assert lammps_task.structure_filename == structure_filename
    assert isinstance(lammps_task.structure,crystal.SimulationCell)

    assert lammps_task.status == 'INIT'

def test__configure_potential():
    #<--- expected behavior
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
    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
            task_name = task_name,
            task_directory = task_directory,
            structure_filename = structure_filename)
    lammps_task.configuration = configuration
    #<--- test setup
    assert os.path.exists(task_directory)
    assert lammps_task.status == 'INIT'
    assert lammps_task.potential is None
    assert lammps_task.configuration['potential']['potential_type'] \
            == configuration['potential']['potential_type']
        
    #<--- code being testing
    lammps_task.configure_potential()
   
    #<--- testing to see if the potential is setup correctly
    assert isinstance(lammps_task.potential,potential.Potential)
    if type(lammps_task.potential) == potential.EamPotential:
        assert isinstance(
                lammps_task.potential.obj_pair,
                potential.PairPotential)
        assert isinstance(
                lammps_task.potential.obj_density,
                potential.EamDensityFunction)
        assert isinstance(
                lammps_task.potential.obj_embedding,
                potential.EamEmbeddingFunction)

    if type(lammps_task.potential) == potential.EamPotential \
            and 'setfl_filename' in configuration['potential']:
        if configuration['potential']['setfl_filename'] is not None:
            lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))
            assert os.path.isfile(os.path.join(
                    task_directory,
                    lammps_setfl_filename))

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
    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
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

    if len(lammps_task.conditions_READY) == 0:
        assert lammps_task.status == 'READY'
    else:
        assert lammps_task.status == 'CONFIG'

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
    #<--- expected behavior
    lammps_setfl_filename = '{}.eam.alloy'.format("".join(symbols))

    #<--- check to see if the configuration dictionaries were setup correctly
    assert type(lammps_task.conditions_INIT) == OrderedDict
    assert type(lammps_task.conditions_CONFIG) == OrderedDict
    assert type(lammps_task.conditions_READY) == OrderedDict
    assert type(lammps_task.conditions_RUNNING) == OrderedDict
    assert type(lammps_task.conditions_POST) == OrderedDict
    assert type(lammps_task.conditions_FINISHED) == OrderedDict

    assert os.path.isfile(os.path.join(
        task_directory,
        lammps_setfl_filename))

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



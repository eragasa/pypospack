import pytest
import os,shutil,copy
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
        'test_LammpsElasticCalculation',
        'MgO_NaCl_unit.gga.relax.vasp')

MgO_LC_configuration = OrderedDict()
MgO_LC_configuration['task'] = OrderedDict()
MgO_LC_configuration['task']['task_name'] = 'MgO_NaCl.E_sp'
MgO_LC_configuration['task']['task_directory'] = 'MgO_NaCl.E_sp'
MgO_LC_configuration['task_type'] = 'min_none'
MgO_LC_configuration['potential'] = MgO_buck_potential_definition
MgO_LC_configuration['parameters'] = MgO_LC_parameters
MgO_LC_configuration['structure'] = MgO_structure_definition

configuration=MgO_LC_configuration

def cleanup(task_directory):
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsSimulation

def test____init___():
    task_name = configuration['task']['task_name']
    task_directory = configuration['task']['task_directory']
    structure_filename = configuration['structure']['filename']
    restart=False
    fullauto=False

    cleanup(task_directory)
    assert not os.path.exists(task_directory)
    #<--- code being tested
    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
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
    # assert lammps_task.lammps_eam_filename is None
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
    # assert lammps_task.lammps_eam_filename is None
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

import pypospack.task.lammps as tsk_lammps
import pypospack.potential
import pytest

"""
This file contains testing and integration code for elastic calculations
using the buckingham interatomic potential.

to run the unit and integration tests.
>>> pytest

to run sample code
>>> python test_task_LammpsElasticCalculation.py

"""

def get_lewis_catlow_param_dict():
    param_dict = {}
    param_dict['chrg_Mg'] = +2.0
    param_dict['chrg_O']  = -2.0
    param_dict['MgMg_A']   = 0.0 
    param_dict['MgMg_rho'] = 0.5
    param_dict['MgMg_C']   = 0.0
    param_dict['MgO_A']    = 821.6
    param_dict['MgO_rho']  = 0.3242
    param_dict['MgO_C']    = 0.0
    param_dict['OO_A']     = 2274.00 
    param_dict['OO_rho']   = 0.1490
    param_dict['OO_C']     = 27.88

    return copy.deepcopy(param_dict)

class TestTaskLammpsElasticCalculation(object):

    @classmethod 
    def setup_class(self):
        self.task_name = 'task_name'
        self.task_directory = 'task_directory'

        self.structure_db = 'test_LammpsElasticCalculation'
        self.structure_filename = 'MgO_NaCl_unit.vasp'

        self.structure_dict = {}
        self.structure_dict['name'] = 'MgO_NaCl_unit'
        self.structure_dict['filename'] = os.path.join(\
                self.structure_db, 
                self.structure_filename)

        self.param_dict = {}
        self.param_dict['chrg_Mg'] = +2.0
        self.param_dict['chrg_O']  = -2.0
        self.param_dict['MgMg_A']   = 0.0 
        self.param_dict['MgMg_rho'] = 0.5
        self.param_dict['MgMg_C']   = 0.0
        self.param_dict['MgO_A']    = 821.6
        self.param_dict['MgO_rho']  = 0.3242
        self.param_dict['MgO_C']    = 0.0
        self.param_dict['OO_A']     = 2274.00 
        self.param_dict['OO_rho']   = 0.1490
        self.param_dict['OO_C']     = 27.88

        self.potential_dict = {}
        self.potential_dict['potential_type'] = 'buckingham'
        self.potential_dict['elements'] = ['Mg','O']
        self.potential_dict['params'] = copy.deepcopy(self.param_dict)

        self.config_dict = {}
        self.config_dict['structure'] = copy.deepcopy(self.structure_dict)
        self.config_dict['potential'] = copy.deepcopy(self.potential_dict)

        if os.path.exists(self.task_directory):
            shutil.rmtree(self.task_directory)

    def class_init(self):
        self.task = tsk_lammps.LammpsElasticCalculation(
                self.task_name, 
                self.task_directory)

    def class_config(self):
        self.task.config(
                structure = self.structure_dict,
                potential = self.potential_dict)

    def class_ready(self):
        self.task.ready()

    def class_run(self):
        self.task.run(is_mpi=False)

    def class_post(self):
        self.task.postprocess()

    def test_init(self):
        self.task = tsk_lammps.LammpsElasticCalculation(
                self.task_name, 
                self.task_directory)
        assert os.path.exists(self.task_directory)
        assert self.task.task_name == self.task_name
        assert self.task.task_directory == os.path.join(\
                os.getcwd(),
                self.task_directory)
        assert self.task.status == 'INIT'

    def test_config(self):
        self.class_init()
        self.class_config()
        assert self.task.status == 'CONFIG'
        assert isinstance(self.task.potential,pypospack.potential.Buckingham)

    def test_ready(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        assert self.task.status == 'READY'

    def test_run(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        self.class_run()
        assert self.task.status == 'POST'

    def test_post(self):
        self.class_init()
        self.class_config()
        self.class_ready()
        self.class_run()
        self.class_post()

if __name__ == "__main__":
    task_name = 'task_name'
    task_directory = 'task_name'

    structure_db = 'rsrc'
    structure_filename = 'MgO_NaCl_unit.vasp'

    structure_dict = {}
    structure_dict['name'] = 'MgO_NaCl_unit'
    structure_dict['filename'] = os.path.join(\
            structure_db, structure_filename)

    config_dict = {}
    config_dict['structure'] = {}
    config_dict['structure']['name'] = 'MgO_NaCl_unit'
    config_dict['structure']['filename'] = 'rsrc/MgO_NaCl_unit.vasp'
    config_dict['filename'] = structure_filename
    config_dict['potential'] = {}
    config_dict['potential']['potential_type'] = 'buckingham'
    config_dict['potential']['elements'] = ['Mg','O']
    config_dict['potential']['params'] = get_lewis_catlow_param_dict()

    print(80*'-')
    print('task.status: null - > INIT')
    task = tsk_lammps.LammpsElasticCalculation(
            task_name = task_name,
            task_directory = task_directory)

    print(80*'-')
    print('task.status: INIT -> CONFIG')
    if task.status == 'INIT':
        task.config(
                structure=config_dict['structure'],
                potential=config_dict['potential'])

    print(80*'-')
    print('task.status: CONFIG -> READY')
    if task.status == 'CONFIG':
        pre_var_dict = {}
        task.ready(pre_var_dict)

    print(80*'-')
    print('task.status: READY -> RUN')
    if task.status == 'READY':
        param_dict = get_lewis_catlow_param_dict()
        task.run(param_dict)

    print(80*'-')
    print('task.status: READY -> POST')
    if task.status == 'RUN':
        pass

    print(80*'-')
    print('task.status: POST -> DONE')
    if task.status == 'POST':
        task.post()
        print('task.results')
        print(task.results)

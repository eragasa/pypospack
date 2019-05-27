import os, copy, shutil, subprocess
from collections import OrderedDict
import pypospack.crystal as crystal
import pypospack.io.vasp as vasp
import pypospack.potential as potential

def test__import__from_pypospack_task_gulp():
    from pypospack.task.gulp import GulpPhononCalculation

MgO_LewisCatlow= OrderedDict()
MgO_LewisCatlow['chrg_Mg'] = +2.0
MgO_LewisCatlow['chrg_O']  = -2.0
MgO_LewisCatlow['MgMg_A']   = 0.0 
MgO_LewisCatlow['MgMg_rho'] = 0.5
MgO_LewisCatlow['MgMg_C']   = 0.0
MgO_LewisCatlow['MgO_A']    = 821.6
MgO_LewisCatlow['MgO_rho']  = 0.3242
MgO_LewisCatlow['MgO_C']    = 0.0
MgO_LewisCatlow['OO_A']     = 2274.00 
MgO_LewisCatlow['OO_rho']   = 0.1490
MgO_LewisCatlow['OO_C']     = 27.88


configuration_MgO = OrderedDict()
configuration_MgO['potential'] = OrderedDict()
configuration_MgO['potential']['potential_type']='buckingham'
configuration_MgO['potential']['symbols'] = ['Mg','O']
configuration_MgO['potential']['parameter_names']=[
        'chrg_Mg','chrg_O',
        'MgMg_A','MgMg_rho','MgMg_C',
        'MgO_A','MgO_rho','MgO_C',
        'OO_A','OO_rho','OO_C']
configuration_MgO['parameters'] = copy.deepcopy(MgO_LewisCatlow)

def test__init__():
    #<--- variable setup
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False

    #<--- file directory setup
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    #<--- code being tested
    from pypospack.task.gulp import GulpPhononCalculation
    gulp_task = GulpPhononCalculation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)

    #<--- test
    assert gulp_task.task_name == task_name
    assert os.path.abspath(gulp_task.task_directory) \
            == os.path.abspath(task_directory)
    assert gulp_task.root_directory == os.getcwd()
    assert os.path.exists(task_directory)
    assert gulp_task.structure_filename == structure_filename
    assert isinstance(
            gulp_task.structure,
            crystal.SimulationCell)
    assert gulp_task.gulp_input_filename == 'gulp.in'
    assert gulp_task.gulp_output_filename == 'gulp.out'
    assert gulp_task.gulp_bin == os.environ['GULP_BIN']
    assert gulp_task.status == 'INIT'

    #<--- cleanup
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

def test__write_gulp_input_file():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)
    #<--- test of the setup variables
    assert os.path.isfile(structure_filename)
    #<--- setup of the file system
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)
    #<--- code setup
    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)
    gulp_task.configuration = configuration
    gulp_task.read_structure_file()
    gulp_task.configure_potential()

    #<--- test preconditions
    assert not os.path.isfile(os.path.join(
        gulp_task.task_directory,
        gulp_task.gulp_input_filename))

    #<--- code being tested
    gulp_task.write_gulp_input_file()
    
    #<--- expected results
    symbols = configuration['potential']['symbols']
    parameter_names = configuration['potential']['parameter_names']
    
    #<--- testing the results
    assert os.path.isfile(os.path.join(
        gulp_task.task_directory,
        gulp_task.gulp_input_filename))

def test__configure_potential():
    #<--- setup of variables
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)
    #<--- some tests for the variable setup 
    assert 'potential' in configuration
    assert 'potential_type' in configuration['potential']
    #<--- setup directory structure and filesystem for test
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    #<--- code setup
    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)
    #<--- code being tested
    gulp_task.configure_potential(
            configuration=configuration)
    #<--- expected results
    symbols = configuration['potential']['symbols']
    parameter_names = configuration['potential']['parameter_names']
    #<--- tests
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is None
def test__on_init():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)
    #<--- testing the setup
    assert 'parameters' in configuration
    #<--- setup the filesystem
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)
    #<--- code setup
    from pypospack.task.gulp import GulpPhononCalculation
    gulp_task = GulpPhononCalculation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)

    #<--- expected values
    symbols=configuration['potential']['symbols']
    parameter_names=configuration['potential']['parameter_names']
    
    #<--- precondition test
    assert gulp_task.status == 'INIT'

    #<--- code being tested
    gulp_task.on_init(configuration=configuration)
    
    #<------ test attribute.potential
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    #<------ test attribute.symbol
    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is not None
   
    #<----- check that status update is correct
    assert gulp_task.conditions_CONFIG['is_structure_configured'] is True
    assert gulp_task.conditions_CONFIG['is_potential_configured'] is True
    assert all([v for k,v in gulp_task.conditions_INIT.items()]) 
    assert all([v for k,v in gulp_task.conditions_CONFIG.items()]) 
    assert gulp_task.all_conditions_INIT
    assert gulp_task.all_conditions_CONFIG
    assert gulp_task.all_conditions_INIT and gulp_task.all_conditions_CONFIG
    assert gulp_task.status == 'CONFIG'


def test__on_config():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)

    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    #<--- code setup
    from pypospack.task.gulp import GulpPhononCalculation
    gulp_task = GulpPhononCalculation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)
    gulp_task.configuration = copy.deepcopy(configuration)
    gulp_task.on_init()
    
    #<--- precondition test
    assert gulp_task.status == 'CONFIG'

    #<--- code being tested
    gulp_task.on_config()
    
    #<--- expected values
    symbols=configuration['potential']['symbols']
    parameter_names=configuration['potential']['parameter_names']
    
    #<------ test attribute.potential
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    #<------ test attribute.symbol
    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is not None
   
    #<----- check that status update is correct
    assert gulp_task.conditions_READY['input_file_is_written']
    assert gulp_task.all_conditions_INIT
    assert gulp_task.all_conditions_CONFIG
    assert gulp_task.all_conditions_READY
    assert gulp_task.status == 'READY'

def test__on_ready():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpPhononCalculation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    
    configuration=copy.deepcopy(configuration_MgO)

    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    #<--- code setup
    from pypospack.task.gulp import GulpPhononCalculation
    gulp_task = GulpPhononCalculation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)
    gulp_task.configuration = copy.deepcopy(configuration)
    gulp_task.on_init()
    gulp_task.on_config()
    #<--- precondition test
    assert gulp_task.status == 'READY'

    #<--- code being tested
    gulp_task.on_ready()
    
    #<--- expected values
    symbols=configuration['potential']['symbols']
    parameter_names=configuration['potential']['parameter_names']
    
    #<------ test attribute.potential
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    #<------ test attribute.symbol
    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is not None
   
    #<----- check that status update is correct
    assert gulp_task.conditions_READY['input_file_is_written']
    assert gulp_task.all_conditions_INIT
    assert gulp_task.all_conditions_CONFIG
    assert gulp_task.all_conditions_READY
    assert gulp_task.all_conditions_POST
    assert gulp_task.status == 'POST'



if __name__ == '__main__':
    vasp_filename = 'MgO_NaCl_unit.vasp'
    vasp_input_filename = os.path.join(os.getcwd(),'rsrc',vasp_filename)

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

    #### TEST INITIALIZE ####
    task_name = 'gulp_test'
    task_directory = os.path.join(\
            os.getcwd(),
            task_name)
    gulp_input_filename = 'gulp.in'
    #task = GulpPhononCalculation(task_name,task_directory)

    #### TEST POTENTIAL POTENTIAL ####
    #print('----- test that buckingham potential provides the right format -----')
    #task.potential = potential.Buckingham(['Mg','O'])
    #task.param_dict = copy.deepcopy(param_dict)
    #print(task.potential.gulp_potential_section_to_string(param_dict))

    #### TEST IF WE CAN CAN WRITE THE INPUT FILE ####
    #gulp_input_filename = os.path.join(task.task_directory,'gulp.in')
    #poscar_input_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    #task.write_gulp_input_file(\
    #        filename=gulp_input_filename,
    #        poscar=vasp_input_filename)

    #### TEST IF WE CAN RUN THE BUCKINGHAM POTENTIAL ####
    #task.structure_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    #task.run()


    task = GulpPhononCalculation(task_name,task_directory)
    task.potential = potential.Buckingham(['Mg','O'])
    task.param_dict = copy.deepcopy(param_dict)
    task.structure_file = os.path.join(\
            task.task_directory,
            vasp_input_filename)
    task.write_gulp_input_file(\
            filename=os.path.join(task.task_directory,gulp_input_filename),
            poscar=vasp_input_filename)
    task.run()

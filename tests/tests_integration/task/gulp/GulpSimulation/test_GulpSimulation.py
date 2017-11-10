import pytest 
import os, shutil, copy
from collections import OrderedDict
import pypospack.potential as potential

def test__import__from_pypospack_task_gulp():
    from pypospack.task.gulp import GulpSimulation

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

def test____init__():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpSimulation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False

    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)

    shutil.rmtree(task_directory)

def test__configure_potential():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpSimulation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)

    gulp_task.configure_potential(
            configuration['potential'])

    symbols = configuration['potential']['symbols']
    parameter_names = configuration['potential']['parameter_names']
    
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is None

def test__configure():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpSimulation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    
    configuration=copy.deepcopy(configuration_MgO)

    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    #<--- code setup
    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
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
    gulp_task.configure(configuration=configuration)
    
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
    gulp_task.update_status()
    assert gulp_task.status == 'CONFIG'

def test__write_gulp_input_file():
    task_name = 'gulp_test'
    task_directory = 'gulp_test_directory'
    structure_filename = os.path.join(
            'test_GulpSimulation',
            'MgO_NaCl_prim.gga.relax.vasp')
    restart=False
    configuration=copy.deepcopy(configuration_MgO)
    if os.path.exists(task_directory):
        shutil.rmtree(task_directory)

    from pypospack.task.gulp import GulpSimulation
    gulp_task = GulpSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            restart=restart)

    gulp_task.read_structure_file()
    gulp_task.configure_potential(configuration['potential'])

    symbols = configuration['potential']['symbols']
    parameter_names = configuration['potential']['parameter_names']
    
    assert type(gulp_task.potential) is potential.BuckinghamPotential
    assert gulp_task.potential.symbols == symbols
    assert gulp_task.potential.parameter_names == parameter_names

    assert isinstance(gulp_task.parameters,dict)
    for pn in parameter_names:
        assert pn in gulp_task.parameters
    for pn,pv in gulp_task.parameters.items():
        assert pv is None
if __name__ == '__main__':
    #### TEST INITIALIZE ####
    from pypospack.task.gulp import GulpSimulation
    task = GulpSimulation(
            task_name='gulp_task',
            task_directory='gulp_simulation',
            structure_filename=os.path.join(
                'test_GulpSimulation',
                'MgO_NaCl_prim.gga.relax.vasp'),
            restart=False)

    task.configure(configuration=configuration_MgO)

    #### TEST POTENTIAL POTENTIAL ####
    print('----- test that buckingham potential provides the right format -----')
    task.potential = potential.Buckingham(['Mg','O'])
    task.param_dict = copy.deepcopy(param_dict)
    print(task.potential.gulp_potential_section_to_string(param_dict))

    #### TEST IF WE CAN CAN WRITE THE INPUT FILE ####
    gulp_input_filename = os.path.join(task.task_directory,'gulp.in')
    poscar_input_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    task.write_gulp_input_file(\
            filename=gulp_input_filename,
            poscar=vasp_input_filename)

    #### TEST IF WE CAN RUN THE BUCKINGHAM POTENTIAL ####
    task.structure_filename = os.path.join('rsrc','MgO_NaCl_prim.vasp')
    task.run()


    #task = GulpPhononCalculation(task_name,task_directory)
    #task.potential = potential.Buckingham(['Mg','O'])
    #task.param_dict = copy.deepcopy(param_dict)
    #task.write_gulp_input_file(\
    #        filename=gulp_input_filename,
    #        poscar=vasp_input_filename)
    #task.run()

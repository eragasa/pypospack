import pytest
import os,shutil

def test__import_from_pypospack_task_phonts():
    from pypospack.task.phonts import PhontsSimulation

def test____init____norestart__with_cleanup():
    #<--- create simulation variables
    task_name = "task_name"
    task_directory = "task_directory"
    restart = False
    structure_filename='test_PhontsSimulation/Si_dia_unit.gga.relaxed.vasp'
    phonts_filename='phonons_input.dat'
    root_directory = os.getcwd()

    #<--- setup conditions for the test
    #remove any previous directory that have already existed
    #we are testing for the case where we assume that no
    #directory exists.
    if os.path.isdir(task_directory):
        shutil.rmtree(task_directory)

    #<--- code we are testing
    from pypospack.task.phonts import PhontsSimulation
    testtask = PhontsSimulation(
            task_name=task_name,
            task_directory=task_directory,
            structure_filename=structure_filename,
            phonts_filename=phonts_filename,
            restart=False)

    #<---- testing attributes which should be should be set correctly
    #      from the parent class
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(testtask.task_directory)
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(testtask.task_directory)
    assert testtask.task_name == task_name
    assert os.path.isdir(testtask.task_directory)
    #<---- checking attributes
    assert type(testtask.structure_filename) is str
    assert testtask.structure_filename == structure_filename
    assert os.path.abspath(testtask.phonts_filename) \
            == os.path.abspath(os.path.join(task_directory,phonts_filename))
    assert testtask.structure is None
    assert testtask.potential is None
    assert testtask.interface is None
    assert testtask.fp_interface is None
    assert testtask.delta == 0.005

    #<---- cleanup
    shutil.rmtree(task_directory)

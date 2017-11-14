import pytest
import os,shutil

def test____import_from_pypospack_task_base():
    from pypospack.task.base import Task

def test____import_from_pypospack_task_Task():
    from pypospack.task import Task


def test____init____norestart__with_cleanup():
    #<--- create simulation variables
    task_name = "task_name"
    task_directory = "task_directory"
    restart = False
    root_directory = os.getcwd()
    #<--- setup conditions for the test
    if os.path.isdir(task_directory):
         shutil.rmtree(task_directory)

    #<--- code we are testing
    from pypospack.task import Task
    testtask = Task(
            task_name=task_name,
            task_directory=task_directory,
            restart=False)

    #<--- did the code leave us where we started?
    assert os.getcwd() == root_directory

    #<--- checking the attributes
    assert type(testtask.is_restart) is bool
    assert testtask.is_restart == restart

    assert type(testtask.root_directory) is str
    assert testtask.root_directory == os.getcwd()

    assert type(testtask.task_directory) is str
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(task_directory)
    assert testtask.task_name == task_name

    #<--- check directory structure
    assert os.path.isdir(testtask.task_directory)
    #<--- cleanup
    shutil.rmtree(task_directory)


def test____init____norestart__no_cleanup():
    #<--- create simulation variables
    task_name = "task_name"
    task_directory = "task_directory"
    restart = False
    root_directory = os.getcwd()
    #<--- setup conditions for the test
    if os.path.isdir(task_directory):
         shutil.rmtree(task_directory)
    #<------- this creates an existing simulation directory
    from pypospack.task import Task
    testtask = Task(
            task_name=task_name,
            task_directory=task_directory,
            restart=False)
    assert os.path.isdir(task_directory)
    #<--- code we are testing
    try:
        testtask = Task(
                task_name=task_name,
                task_directory=task_directory,
                restart=False)
    except:
        pytest.fail()
    #<--- did the code leave us where we started?
    assert os.getcwd() == root_directory

    #<--- checking the attributes
    assert type(testtask.is_restart) is bool
    assert testtask.is_restart == restart

    assert type(testtask.root_directory) is str
    assert testtask.root_directory == os.getcwd()

    assert type(testtask.task_directory) is str
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(task_directory)
    assert testtask.task_name == task_name

    #<-- check directory structure
    assert os.path.isdir(testtask.task_directory)

    #<--- cleanup
    shutil.rmtree(task_directory)

def test____init____restart__with_cleanup():
    #<--- create simulation variables
    task_name = "task_name"
    task_directory = "task_directory"
    restart = True
    root_directory = os.getcwd()
    #<--- setup conditions for the test
    if os.path.isdir(task_directory):
         shutil.rmtree(task_directory)

    #<--- code we are testing
    from pypospack.task import Task
    testtask = Task(
            task_name=task_name,
            task_directory=task_directory,
            restart=restart)

    #<--- did the code leave us where we started?
    assert os.getcwd() == root_directory

    #<--- checking the attributes
    assert testtask.is_restart == restart
    assert testtask.root_directory == os.getcwd()
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(task_directory)
    assert testtask.task_name == task_name

    #<-- check directory structure
    assert os.path.isdir(testtask.task_directory)
    #<--- cleanup
    shutil.rmtree(task_directory)

def test__init____restart__no_cleanup():
    #<--- create simulation variables
    task_name = "task_name"
    task_directory = "task_directory"
    restart = True
    root_directory = os.getcwd()
    #<--- setup conditions for the test
    if os.path.isdir(task_directory):
         shutil.rmtree(task_directory)
    #<------- this creates an existing simulation directory
    from pypospack.task import Task
    testtask = Task(
            task_name=task_name,
            task_directory=task_directory,
            restart=restart)
    assert os.path.isdir(task_directory)
    #<--- code we are testing
    testtask = Task(
            task_name=task_name,
            task_directory=task_directory,
            restart=restart)
    #<--- did the code leave us where we started?
    assert os.getcwd() == root_directory
    #<--- checking the attributes
    assert testtask.is_restart == restart
    assert testtask.root_directory == os.getcwd()
    assert os.path.abspath(testtask.task_directory) \
            == os.path.abspath(task_directory)
    assert testtask.task_name == task_name
    #<-- check directory structure
    assert os.path.isdir(testtask.task_directory)
    #<--- cleanup
    shutil.rmtree(task_directory)


import pytest
from collections import OrderedDict
from pypospack.qoi import QoiDatabase
from pypospack.qoi import QoiManager
from pypospack.task import TaskManager
from pypospack.task import Task
from pypospack.pyposmat import PyposmatDataFile
from pypospack.pyposmat import PyposmatEngine
from pypospack.pyposmat import PyposmatConfigurationFile
from MgO import MgO_LewisCatlow

testid1 = 'MgO_NaCl.Ecoh'
testcase1 = ('MgO_NaCl.Ecoh','Ecoh_min_all',OrderedDict([('ideal','MgO_NaCl')]),4.5)

testlabels = "qoi_name,qoi_type,structures,target"
testdata = [testcase1]
testids = [testid1]


def write_configuration_file(qoi_name,qoi_type,structures,target):
    qoidb = QoiDatabase()
    qoidb.add_qoi(
            qoi_name=qoi_name,
            qoi_type=qoi_type,
            structures=structures,
            target=target
        )
    
    potential = OrderedDict()
    potential['potential_type'] = 'buckingham'
    potential['symbols'] = ['Mg','O']
    potential['cutoff_global'] = 10.0
    
    structures = OrderedDict()
    structures['structure_directory'] = 'test_PypospackEngine'
    structures['structures'] = OrderedDict()
    structures['structures']['MgO_NaCl'] = 'MgO_NaCl_unit.gga.relax.vasp'
    
    configuration = PyposmatConfigurationFile()
    configuration.qois = qoidb.qois
    configuration.potential = potential
    configuration.structures = structures
    configuration.write(filename='pyposmat.config.in')

@pytest.mark.parametrize("qoi_name,qoi_type,structures,target",testdata,ids=testids)
def test__init__w_filename(qoi_name,qoi_type,structures,target):
    _filename='pyposmat.config.in'
    write_configuration_file(qoi_name,qoi_type,structures,target)

    #<------------- code being tested
    pyposmatengine = PyposmatEngine(filename_in=_filename,fullauto=False)

    # <------------ test for expected results
    assert pyposmatengine.pyposmat_filename_in == _filename
    assert isinstance(pyposmatengine.configuration,type(None))
    assert isinstance(pyposmatengine.task_manager,type(None))
    assert isinstance(pyposmatengine.qoi_manager,type(None))

@pytest.mark.parametrize("qoi_name,qoi_type,structures,target",testdata,ids=testids)
def test__read_configuration_file(qoi_name,qoi_type,structures,target):
    # <------------ setup for this test
    _filename='pyposmat.config.in'
    write_configuration_file(qoi_name,qoi_type,structures,target)
    pyposmatengine = PyposmatEngine(filename_in=_filename,fullauto=False)
    
    # <------------ code being tested
    pyposmatengine.read_configuration_file()

    # <------------ test for expected results
    assert isinstance(pyposmatengine.configuration,PyposmatConfigurationFile)
    assert isinstance(pyposmatengine.task_manager,type(None))
    assert isinstance(pyposmatengine.qoi_manager,type(None))

@pytest.mark.parametrize("qoi_name,qoi_type,structures,target",testdata,ids=testids)
def test__configure_qoi_manager(qoi_name,qoi_type,structures,target):
    # <------------ setup for this test
    _filename='pyposmat.config.in'
    write_configuration_file(qoi_name,qoi_type,structures,target)
    pyposmatengine = PyposmatEngine(filename_in=_filename,fullauto=False)
    pyposmatengine.read_configuration_file()
    
    # <------------ code being tested
    pyposmatengine.configure_qoi_manager()

    # <------------ testing for expected results
    assert isinstance(pyposmatengine.configuration,PyposmatConfigurationFile)
    assert isinstance(pyposmatengine.qoi_manager,QoiManager)
    assert isinstance(pyposmatengine.task_manager,type(None))

@pytest.mark.parametrize("qoi_name,qoi_type,structures,target",testdata,ids=testids)
def test__configure_task_manager(qoi_name,qoi_type,structures,target):
    # <------------ setup for this test
    _filename='pyposmat.config.in'
    write_configuration_file(qoi_name,qoi_type,structures,target)
    pyposmatengine = PyposmatEngine(filename_in=_filename,fullauto=False)
    pyposmatengine.read_configuration_file()
    pyposmatengine.configure_qoi_manager()

    # <------------ code being tested
    pyposmatengine.configure_task_manager()

    # <------------ testing for expected results
    assert isinstance(pyposmatengine.configuration,PyposmatConfigurationFile)
    assert isinstance(pyposmatengine.qoi_manager,QoiManager)
    assert isinstance(pyposmatengine.task_manager,TaskManager)
    assert isinstance(pyposmatengine.task_manager.tasks,OrderedDict)
    assert isinstance(pyposmatengine.task_manager.obj_Task,OrderedDict)
    for k_task,v_task in pyposmatengine.task_manager.obj_Task.items():
        print(k_task,':',type(v_task),':',v_task.status)
    assert all([isinstance(v_task,Task) 
        for k_task,v_task in pyposmatengine.task_manager.obj_Task.items()])
    assert all([v_task.status == 'INIT'
        for k_task,v_task in pyposmatengine.task_manager.obj_Task.items()])

@pytest.mark.parametrize("qoi_name,qoi_type,structures,target",testdata,ids=testids)
def test__configure_task_manager(qoi_name,qoi_type,structures,target):
    # <------------ setup for this test
    _filename='pyposmat.config.in'
    write_configuration_file(qoi_name,qoi_type,structures,target)
    pyposmatengine = PyposmatEngine(filename_in=_filename,fullauto=False)
    pyposmatengine.read_configuration_file()
    pyposmatengine.configure_qoi_manager()
    pyposmatengine.configure_task_manager()

    # <------------ code being tested
    pyposmatengine.evaluate_parameter_set(
            parameters=MgO_LewisCatlow['parameters'])


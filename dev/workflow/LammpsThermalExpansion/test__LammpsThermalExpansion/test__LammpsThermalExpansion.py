import pytest
import os, shutil
from collections import OrderedDict

import pypospack.utils

from pypospack.task.lammps import LammpsNptSimulation
from pypospack.workflows import LammpsThermalExpansion
# Stillinger and Weber, Phy. Rev. B. v. 31, p.5262 (1985)
Si__stillingerweber = OrderedDict()
Si__stillingerweber['potential_definition'] = OrderedDict()
Si__stillingerweber['potential_definition']['potential_type'] = 'stillingerweber'
Si__stillingerweber['potential_definition']['symbols'] = ['Si']
Si__stillingerweber['parameters'] = OrderedDict([
    ('SiSiSi_epsilon',2.1686),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',21.0),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])

# the structures required for the simulation, which would come from the simulation database
Si_structure_definition = OrderedDict()
Si_structure_definition['name'] = 'Si_dia_unit'
Si_structure_definition['filename'] = os.path.join(
        os.path.join(
            pypospack.utils.get_pypospack_root_directory(),
            'data','Si__structure_db',
            'Si_dia_unit.vasp')
)

workflow_definition = OrderedDict()
workflow_definition['temperature_low'] = 100
workflow_definition['temperature_high'] = 1000
workflow_definition['temperature_step'] = 100
workflow_definition['pressure'] = 0
workflow_definition['time_step'] = 0.001
workflow_definition['time_total'] = 0.001 * 10000
workflow_definition['supercell'] = [10,10,10]

def dev____init__():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    print(workflow.workflow_type)
    print(workflow.workflow_name)
    print(workflow.workflow_path)
    print(workflow.temperatures)

def dev__create_task_configurations():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    for task_name,task_config in workflow.task_configurations.items():
        print(task_name)
        for k,v in task_config.items():
            print('\t',k,v)

def teardown_function():
    path_to_delete = 'Si_dia_unit.thermal_expansion'
    
    if os.path.isdir(path_to_delete):
        shutil.rmtree(path_to_delete)

def test__create_tasks():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    workflow.create_tasks()

def test____init__():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    assert workflow.workflow_type == LammpsThermalExpansion.workflow_type
    assert workflow.workflow_name == LammpsThermalExpansion.get_workflow_name(
            structure_name = Si_structure_definition['name'])
    assert os.path.isdir(workflow.workflow_path)
    assert isinstance(workflow.temperatures,list)
    assert len(workflow.temperatures) == 10
    assert workflow.temperatures == [100,200,300,400,500,600,700,800,900,1000]

def test__create_task_configurations():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    for task_name,task_config in workflow.task_configurations.items():
        print(task_name)
        for k,v in task_config.items():
            print('\t',k,v)

def test__create_tasks():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    workflow.create_tasks()

    for task_name,task_object in workflow.tasks.items():
        assert isinstance(task_object,LammpsNptSimulation)

def test__prepare_tasks():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    workflow.create_tasks()
    workflow.prepare_tasks(
            potential_definition = Si__stillingerweber['potential_definition'],
            potential_parameters = Si__stillingerweber['parameters'])

@pytest.mark.long
def test__run():
    workflow = LammpsThermalExpansion(
            structure_name=Si_structure_definition['name'],
            structure_path=Si_structure_definition['filename'],
            **workflow_definition)
    workflow.create_task_configurations()
    workflow.create_tasks()
    workflow.prepare_tasks(
            potential_definition = Si__stillingerweber['potential_definition'],
            potential_parameters = Si__stillingerweber['parameters'])
    workflow.run()


import pytest
import os
from collections import OrderedDict

def test__import__from_pypospack_task_lammps():
    from pypospack.task.lammps import LammpsSimulation


potential_definition = OrderedDict()
potential_definition['potential_type'] = 'buckingham'
potential_definition['symbols'] = ['Mg','O']

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

structure_definition = OrderedDict()
structure_definition['name'] = 'MgO_NaCl_unit'
structure_definition['filename'] = os.path.join(
        'test_LammpsSimulation','MgO_NaCl_unit.gga.relax.vasp')

configuration_MgO = OrderedDict()
configuration_MgO['potential'] = potential_definition
configuration_MgO['parameters'] = MgO_LC_parameters
configuration_MgO['structure'] = structure_definition

def test____init___():
    task_name = 'test_task_name'
    task_directory = 'test_task_directory'

    from pypospack.task.lammps import LammpsSimulation
    lammps_task = LammpsSimulation(
            task_name = task_name,
            task_directory = task_directory)

    #<--- check directory structure
    assert os.path.exists(task_directory)
    assert os.path.abspath(lammps_task.task_directory)\
            == os.path.abspath(task_directory)
    assert lammps_task.status == 'INIT'

    lammps_task.on_init(configuration=configuration_MgO)
    assert lammps_task.status == 'CONFIG'
    assert isinstance(lammps_task.potential,
            pypospack.potential.Potential)

    lammps_task.on_config(configuration=configuration_MgO)
    assert lammps_task.status == 'READY'

    lammps_task.on_ready(configuration=configuration_MgO)
    assert lammps_task.status == 'RUNNING'

    lammps_task.on_running(configuration=configuration_MgO)
    assert lammps_task.status == 'POST'

    lammps_task.on_post(configuration=configuration_MgO)
    assert lammps_task.status == 'FINISHED'



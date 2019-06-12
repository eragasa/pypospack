import pytest
from pypospack.task.lammps import LammpsNptSimulation
pypospack.task.tasks_lammps.lammps_npt_simulation import LammpsNptSimulation


def test__get_task_name():

    task_name = LammpsNptSimulation.get_task_name(structure='structure',temperature=10)
    assert task_name == 'structure.lmps_npt_10'


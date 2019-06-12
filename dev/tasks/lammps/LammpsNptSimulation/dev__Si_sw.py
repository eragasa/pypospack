import sys
import os
import time
import importlib
from collections import OrderedDict

import pypospack.utils
from pypospack.task.lammps import LammpsNptSimulation

# definition of the potential
reference_potentials = OrderedDict()
# Stillinger and Weber,  Phys. Rev. B, v. 31, p. 5262, (1985)
reference_potentials['SW'] = OrderedDict([
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
reference_potentials['VBWM'] = OrderedDict([
    ('SiSiSi_epsilon',1.64833),
    ('SiSiSi_sigma',2.0951),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.5),
    ('SiSiSi_gamma',1.20),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',7.049556277),
    ('SiSiSi_B',0.6022245584),
    ('SiSiSi_p',4.0),
    ('SiSiSi_q',0.0),
    ('SiSiSi_tol',0.0)
])
reference_potentials['PG'] = OrderedDict([
    ('SiSiSi_epsilon',1.04190),
    ('SiSiSi_sigma',2.128117),
    ('SiSiSi_a',1.80),
    ('SiSiSi_lambda',31.0),
    ('SiSiSi_gamma',1.10),
    ('SiSiSi_costheta0',-1/3),
    ('SiSiSi_A',19.0),
    ('SiSiSi_B',0.65),
    ('SiSiSi_p',3.5),
    ('SiSiSi_q',0.5),
    ('SiSiSi_tol',0.0)
])

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


def get_task_name(structure_name,temperature,pressure=0):
    PYPOSPACK_TASK_FMT = "{}.npt_T{}_P{}"

    T_ = int(temperature)
    P_ = int(pressure)
    return PYPOSPACK_TASK_FMT.format(structure_name,T_,P_)

def get_task_directory(structure_name,temperature,pressure=0,base_path="./"):
    task_name_ = get_task_name(structure_name=structure_name, temperature=temperature, pressure=pressure)

    return os.path.join(base_path,task_name_)

def get_npt_task_configuration(
        structure_name, structure_path, temperature, pressure, time_step, time_total,
        supercell=[10,10,10], base_path="./"):
    """
    Args:
       structure_name (str): alias for the structure
       structure_filename (str): path for the structure
       temperature (int): temperature for the simulation in Kelvin
       pressure (int): the pressure for the simualtion in bars
       time_step (int): the time step for the time integrator in picoseconds
       time_total(int): the total time for time integration is picoseconds
       supercell (list): the number of replications in the a1, a2, and a3 dimensions
       base_path (list): the path in which to create the simulations
    """

    assert isinstance(structure_name, str)
    assert os.path.isfile(structure_path)
    assert isinstance(supercell, list)
    assert all([isinstance(v, int) for v in supercell])
    assert len(supercell) == 3

    npt_task_configuration = OrderedDict()
    npt_task_configuration['task_name'] = get_task_name(structure_name,temperature,pressure)
    npt_task_configuration['task_directory'] = get_task_directory(structure_name,temperature,pressure,base_path)
    npt_task_configuration['structure_filename'] = structure_path
    npt_task_configuration['supercell'] = supercell
    npt_task_configuration['temperature'] = temperature
    npt_task_configuration['pressure'] = pressure
    npt_task_configuration['time_step'] = time_step
    npt_task_configuration['time_total'] = time_total
    npt_task_configuration['supercell'] = supercell

    return npt_task_configuration


import inspect
class LammpsNptSimulationScaffold(object):
    """ a development scafford for LammpsNptSimulation

    Tasks are usually run within the TaskManager.  This scafford provides a standard
    interface for running tasks

    """
    def __init__(self,
                 potential_definition,
                 structures_definition,
                 task_configuration,
                 parameters):

        self.potential_definition = potential_definition 
        self.structure_definition = structures_definition
        self.task_configuration = task_configuration
        self.parameters = parameters

        self.lammps_task = LammpsNptSimulation(**task_configuration)

    def run(self,task_configuration=None): 

        if task_configuration is not None:
            self.task_configuration = task_configuration
        self.on_init(self.task_configuration)
        self.on_config(self.task_configuration)
        self.on_ready(self.task_configuration)

    def on_init(self,task_configuration):
        self.task_configuration['potential'] = self.potential_definition
        self.task_configuration['parameters'] = self.parameters
        self.lammps_task.on_init(task_configuration)

    def on_config(self,task_configuration):
        self.lammps_task.on_config(task_configuration)

    def on_ready(self,task_configuration):
        self.lammps_task.on_ready(task_configuration)

        time_sleep = 0.1
        while (self.lammps_task.process.poll() is None):
            time.sleep(time_sleep)
        lammps_task = None

if __name__ == "__main__":
    npt_task_configuration = get_npt_task_configuration(
            structure_name = Si_structure_definition['name'],
            structure_path = Si_structure_definition['filename'],
            time_step = 0.001,
            time_total = 0.001 * 1000,
            pressure = 0,
            temperature = 1000,
            supercell = [10,10,10])

    o = LammpsNptSimulationScaffold(
            potential_definition = Si__stillingerweber['potential_definition'],
            structures_definition = Si_structure_definition,
            task_configuration = npt_task_configuration,
            parameters = Si__stillingerweber['parameters'])
    o.on_init(o.task_configuration)
    o.on_config(o.task_configuration)
    o.on_ready(o.task_configuration)

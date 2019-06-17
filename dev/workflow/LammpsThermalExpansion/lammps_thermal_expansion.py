import sys
import os, shutil
import time
import importlib
from collections import OrderedDict
import numpy as np
import pypospack.utils
from pypospack.task.lammps import LammpsNptSimulation

class Workflow():

    def get_workflow_name(structure_name):
        msg = "subclasses of workflow need to implement this method"
        raise NotImplementedError(msg)

    def create_tasks(self):
        msg = "subclasses of workflow need to implement this method"
        raise NotImplementedError(msg)

    def create_workflow_path(self,path):
        if os.path.isdir(path):
            assert path is not os.getcwd()
            shutil.rmtree(path)
        os.mkdir(path)

class LammpsThermalExpansion(Workflow):
    workflow_type = 'lmps_thermal_expansion'
    def get_workflow_name(structure_name):
        assert isinstance(structure_name,str)
        WORKFLOW_NAME_FORMAT = '{}.thermal_expansion'

        return WORKFLOW_NAME_FORMAT.format(structure_name)
    def __init__(self,
                 structure_name,
                 structure_path,
                 temperature_low,
                 temperature_high,
                 temperature_step,
                 pressure = 0,
                 time_step = 0.001,
                 time_total = 0.001 * 10000,
                 supercell=[10,10,10],
                 workflow_path=None,
                 is_new_workflow=True):

        assert isinstance(structure_name,str)
        self.structure_name = structure_name
        self.structure_path = structure_path
        self.potential_definition = None

        assert isinstance(temperature_low,int) or isinstance(temperature,float)
        assert isinstance(temperature_high,int) or isinstance(temperatuer,float)
        assert isinstance(temperature_step,int) or isinstance(temperature,float)
        assert temperature_low > 0
        assert temperature_high > 0
        assert temperature_step > 0
        assert temperature_low < temperature_high
        self.temperature_low = temperature_low
        self.temperature_high = temperature_high
        self.temperature_step = temperature_step
        self.temperatures = self.get_temperatures()

        assert isinstance(pressure,int) or isinstance(pressure,float)
        assert pressure >= 0
        self.pressure = pressure

        assert isinstance(time_step,int) or isinstance(time_step,float)
        assert isinstance(time_total,int) or isinstance(time_total,float)
        assert time_step > 0
        assert time_total > 0
        self.time_step = time_step
        self.time_total = time_total
        
        assert isinstance(supercell,list)
        assert len(supercell) == 3
        assert all([isinstance(k,int) for k in supercell])
        assert all([k > 0 for k in supercell])
        self.supercell = supercell

        self.workflow_name = LammpsThermalExpansion.get_workflow_name(
                structure_name=self.structure_name)

        if workflow_path is None:
            workflow_path_ = os.path.join(".",self.workflow_name)
        else:
            workflow_path_ = workflow_path
        self.workflow_path = os.path.abspath(workflow_path_)
        if is_new_workflow:
            self.create_workflow_path(path=self.workflow_path)
        else:
            assert os.path.isdir(workflow_path)

        self.tasks = None
        self.task_configurations = None

    def create_task_configurations(self):
        self.task_configurations = OrderedDict()

        for temperature in self.temperatures:
            task_name = self.get_npt_task_name(
                    structure_name=self.structure_name,
                    temperature=temperature,
                    pressure=self.pressure)
            
            task_directory = os.path.join(
                    self.workflow_path,
                    task_name)

            self.task_configurations[task_name] = \
                    LammpsThermalExpansion.get_npt_task_configuration(
                        task_name=task_name,
                        task_directory=task_directory,
                        structure_path=self.structure_path,
                        temperature=temperature,
                        pressure=self.pressure,
                        time_step=self.time_step,
                        time_total=self.time_total,
                        supercell=self.supercell)

    def create_tasks(self):
        self.tasks = OrderedDict()
        for k,v in self.task_configurations.items():
            self.tasks[k] = LammpsNptSimulation(**v)

    def prepare_tasks(self,
                      potential_definition,
                      potential_parameters):
        for task_name,task in self.tasks.items():
            configuration = self.task_configurations[task_name]
            configuration['potential'] = potential_definition
            configuration['parameters'] = potential_parameters
            task.on_init(configuration)
            task.on_config(configuration)

    def run(self):
        for task_name,task in self.tasks.items():
            print(task_name)
            task.on_ready(self.task_configurations[task_name])

            time_sleep = 0.1
            while (task.process.poll() is None):
                time.sleep(time_sleep)

    def prepare_slurm_jobs(self):
        pass


    def get_temperatures(self):

        # determine temperature ranges
        temperatures = []
        temperature = self.temperature_low
        while True:
            temperatures.append(temperature)
            temperature += self.temperature_step
            if temperature > self.temperature_high:
                break

        return temperatures
    def get_npt_task_name(self,structure_name,temperature,pressure):
        PYPOSPACK_TASK_FMT = "{}.npt_T{}_P{}"

        T_ = int(temperature)
        P_ = int(pressure)
        return PYPOSPACK_TASK_FMT.format(structure_name,T_,P_)

    def get_npt_task_configuration(
            task_name, task_directory, structure_path,
            temperature, pressure,
            time_step, time_total,
            supercell=[10,10,10]):
        """
        Args:
            task_name (str): the task name
            task_directory (str): the task directory
            temperature (int,float): the temperature
        """

        assert isinstance(supercell, list)
        assert all([isinstance(v, int) for v in supercell])
        assert len(supercell) == 3

        npt_task_configuration = OrderedDict()
        npt_task_configuration['task_name'] = task_name
        npt_task_configuration['task_directory'] = task_directory
        npt_task_configuration['structure_filename'] = structure_path
        npt_task_configuration['temperature'] = temperature
        npt_task_configuration['pressure'] = pressure
        npt_task_configuration['time_step'] = time_step
        npt_task_configuration['time_total'] = time_total
        npt_task_configuration['supercell'] = supercell

        return npt_task_configuration

    def run_tasks(self):
        pass
    
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
